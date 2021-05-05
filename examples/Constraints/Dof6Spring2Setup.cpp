#include <math.h>
#include <limits>
#include <deque>
#include <tuple>
#include <mutex>
#include <set>

#include "Dof6Spring2Setup.h"

#include "../Utils/b3ResourcePath.h"
#include "Bullet3Common/b3FileUtils.h"
#include "../Importers/ImportObjDemo/LoadMeshFromObj.h"
#include "../OpenGLWindow/GLInstanceGraphicsShape.h"
#include "../Utils/b3BulletDefaultFileIO.h"

#include "BulletDynamics/ConstraintSolver/btGeneric6DofSpring2Constraint.h"
#include "BulletDynamics/Featherstone/btMultiBodyFixedConstraint.h"
#include "BulletDynamics/Featherstone/btMultiBodySphericalJointMotor.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.785398163397448309616
#endif

#include "../CommonInterfaces/CommonRigidBodyBase.h"
#include "../CommonInterfaces/CommonMultiBodyBase.h"

#define LENGTH_RATIO 1

#define MASS_RATIO 10

#define DELTATIME 0.0333333

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////  helper function
////
/////////////////////////////////////////////////////////////
btScalar inline rad2degree(btScalar rad){
    return rad / SIMD_PI * 180.0f;
}

btScalar inline degree2rad(btScalar degree){
    return degree / 180.0f * SIMD_PI;
}

void getEulerValFromQuaternion(btVector3& eulerVal, const btScalar* q)
{
    btQuaternion qt(q[0], q[1], q[2], q[3]);
    qt.getEulerZYX(eulerVal[2], eulerVal[1], eulerVal[0]);
}

btScalar cosOffset(btScalar amp, btScalar T, btScalar phase, btScalar time){
    return amp * cos(time / T * SIMD_2_PI + phase);
}

btScalar sinOffset(btScalar amp, btScalar T, btScalar phase, btScalar time){
    return amp * sin(time / T * SIMD_2_PI + phase);
}

// TODO use it in code
btScalar adjustKs(btScalar ks, btScalar mass, btScalar deltaTime)
{
    btScalar r = ks;
    btScalar angular_freq = btSqrt(r / mass);
    if (0.25 < angular_freq * deltaTime)
    {
        r = BT_ONE / deltaTime / deltaTime / btScalar(16.0) * mass;
    }
    return r;
}

// TODO adjust kd


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////  multibody skeleton
////
/////////////////////////////////////////////////////////////
#define WIRE_FRAME 0

#define POS_ARRAY_SIZE 3

#define ACC_ARRAY_SIZE 1

char mfileName[1024];

const btVector3 zeroTransOffset(0.0, 0.0, 0.0);

enum CollisionTypes
{
    NOTHING  = 0,        //< Collide with nothing
    BONE_BODY = 1 << 1,
    TARGET_BODY = 1 << 2
};

static btScalar radius(0.05);

#define MAX_DRAG_FORCE SIMD_INFINITY
#define MAX_SPRING_FORCE 0.5f

class GlobalState
{
public:
    GlobalState(GlobalState& other) = delete;

    void operator=(const GlobalState&) = delete;

    static GlobalState* GetInstance()
    {
        if (_instance == nullptr) {
            _instance = new GlobalState();
        }
        return _instance;
    }

    void insertCollisionEvent(int link)
    {
        std::lock_guard<std::mutex> lk(m_mtx);
        collisionLinks.insert(link);
    }

    std::set<int> getCollisionEvent()
    {
        std::lock_guard<std::mutex> lk(m_mtx);
        return collisionLinks;
    }

    void clearCollisionEvent()
    {
        std::lock_guard<std::mutex> lk(m_mtx);
        collisionLinks.clear();
    }

private:
    GlobalState() {}
    static GlobalState* _instance;
    std::mutex m_mtx;
    std::set<int> collisionLinks;
};

GlobalState* GlobalState::_instance = nullptr;

ATTRIBUTE_ALIGNED16(class)
customMultiBody : public btMultiBody
{
public:
    BT_DECLARE_ALIGNED_ALLOCATOR();

    //
    // initialization
    //

    customMultiBody(int n_links,               // NOT including the base
                    btScalar mass,             // mass of base
                    const btVector3& inertia,  // inertia of base, in base frame; assumed diagonal
                    bool fixedBase,            // whether the base is fixed (true) or can move (false)
                    bool canSleep, bool deprecatedMultiDof = true)
            : btMultiBody(n_links, mass, inertia, fixedBase, canSleep, deprecatedMultiDof)
    {
        m_maxOmega = 100.f;
        m_maxOmegaX = m_maxOmega;
    }


    virtual ~customMultiBody() {}

    void setMaxOmegaX(btScalar max) { m_maxOmegaX = max; }
    void setMaxOmega(btScalar max) { m_maxOmega = m_maxOmegaX = max; }

    virtual void applyDeltaVeeMultiDof(const btScalar* delta_vee, btScalar multiplier)
    {
        for (int dof = 0; dof < 6 + getNumDofs(); ++dof)
        {
            m_realBuf[dof] += delta_vee[dof] * multiplier;
            if (dof >= 6 && dof % 3 == 0)
            {
                if ( dof % 3 == 0 )
                {
                    btClamp(m_realBuf[dof], -m_maxOmegaX, m_maxOmegaX);
                }
                else {
                    btClamp(m_realBuf[dof], -m_maxOmega, m_maxOmega);
                }
            }
            else
                btClamp(m_realBuf[dof], -m_maxCoordinateVelocity, m_maxCoordinateVelocity);
        }
    }

protected:
    btScalar m_maxOmegaX;
    btScalar m_maxOmega;
};


class gravityGenerator
{
private:
    int m_status;
    int m_step;
    float m_curr_g;

public:
    float m_gravity;
    float m_tup;
    float m_tdown;

    gravityGenerator(float gravity = -0.1, int tup = 10, int tdown = 10) : m_gravity(gravity),
        m_tup(tup), m_tdown(tdown)
    {
        m_step = 0;
        m_status = 0;
        m_curr_g = 0;
    }

    float getGVal()
    {
        float g = 0.0f;
        if ( m_status  == 1 ){
            g = m_gravity * m_step / m_tup;
            if ( m_step < m_tup )
                m_step += 1;
        } else if ( m_status == -1){
            float g = m_curr_g * ( 1.0 - m_step / m_tdown );
            if ( m_step < m_tdown )
                m_step += 1;
            else{
                m_status = 0;
            }
        }

        if ( m_status != -1 )
            m_curr_g = g;

        return g;
    }

    void startGravity()
    {
        if ( m_status == 0 || m_status == -1) {
            m_step = 0;
        }
        m_status = 1;
    }

    void stopGravity()
    {
        if ( m_status == 0 || m_status == -1 )
            return;
        m_status = -1;
        m_step = 0;
    }
};

struct LineOscilator
{
    LineOscilator()
    {
        amp = 2.0f;
        T = 30.0f;
        phase =  0;
    }

    std::tuple<btVector3, btScalar> getPos(float time, float deltaTime){
        btVector3 pos = btVector3(0, 0, cosOffset(amp, T, phase, time) * exp(-time / (4* T)));
        return std::make_tuple(pos, 0);
    }

    btScalar amp;
    btScalar T;
    btScalar phase;
};


struct CircleOscilator
{
    CircleOscilator()
    {
		T = 60 * DELTATIME;
        omega = SIMD_2_PI / T;
        angle = SIMD_PI / 2.0;
        clockwise = true;
		radius = 0.7523666 * LENGTH_RATIO;
    }

    std::tuple<btVector3, btScalar> getPos(float time, float deltaTime){
        if ( angle > (SIMD_PI / 2.0) )
        {
            clockwise = false;
        }
        if ( angle < (-SIMD_PI / 2.0) )
        {
            clockwise = true;
        }
        btScalar currOmega = omega;
        currOmega = omega * deltaTime * exp(-time / (10 * T));
        angle += (clockwise ? 1 : -1) * currOmega;
//        angle += (clockwise ? 1 : -1) * currOmega * deltaTime;
        return std::make_tuple(btVector3(radius * sin(angle), 0, radius * cos(angle)), currOmega);
    }

    btScalar T;
    btScalar omega;
    btScalar angle;
    bool clockwise;
    btScalar radius;
};

struct Skeleton : public CommonMultiBodyBase
{
    customMultiBody* m_multiBody;
	customMultiBody* m_multiBody1;
	customMultiBody* m_multiBody2;
    btRigidBody* m_collider;
    btRigidBody* m_marker;

    int m_solverType;
    btScalar m_time;
    int m_step;
    int m_numLinks;
	int m_numLinks1;
	int m_numLinks2;
    btAlignedObjectArray<btQuaternion> m_balanceRot;
	btAlignedObjectArray<btQuaternion> m_balanceRot1;
	btAlignedObjectArray<btQuaternion> m_balanceRot2;
    std::vector<btScalar> m_Ks;
    gravityGenerator m_g;
    btScalar m_linearDragEffect;
    btScalar m_centrifugalDragEffect;
    bool m_clockwise;
    btMultiBodyPoint2Point* m_p2p;
    btScalar m_radius;
    CircleOscilator m_move;
    btVector3 m_currVel;
    btVector3 m_currAcc;
    std::deque<btVector3> positions;
    std::deque<btVector3> accs;
    std::vector<btMultiBodySphericalJointMotor*> motors;
    bool use_constraint;
    float m_avgAcc;
    GlobalState* m_state;
    float joint_lengths[6];
	btVector3 linkHalfExtents;

public:
    Skeleton(struct GUIHelperInterface* helper);
    virtual ~Skeleton();
    virtual void initPhysics();
    virtual void stepSimulation(float deltaTime);
    virtual void resetCamera()
    {
		float dist = 0.6 * LENGTH_RATIO;
        float pitch = -20;
        float yaw = -90;
		float targetPos[3] = {0.5, 0, 0.0};
        m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
    }
    void addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const float* joint_lengths);
    void moveCollider(const btVector3& pos);
    void moveMarker(const btVector3& pos);
    void applyBaseLinearDragForce(const btVector3& dir);
    void applyBaseCentrifugalForce(float deltaTime, const btVector3& pos, float omega);
    static void OnInternalTickCallback(btDynamicsWorld* world, btScalar timeStep);
    void getLinearAcc(float deltaTime, int m_step, const btVector3& pos);
    void limitMaxTwist(float max_angle);
	std::pair< btVector3, btVector3> applySpringForce(customMultiBody* mb, const btAlignedObjectArray<btQuaternion>& balanceRot, const btVector3& tailTorque = btVector3(0, 0, 0));
};

Skeleton::Skeleton(struct GUIHelperInterface* helper)
        : CommonMultiBodyBase(helper)
{
    m_time = btScalar(0.0);
    m_step = 1;
    m_numLinks = 6;
	m_numLinks1 = 3;
	m_numLinks2 = 4;
    m_solverType = 0;
    m_g.m_gravity = -0.08;
    m_linearDragEffect = btScalar(25.0);
    m_centrifugalDragEffect = btScalar(0.3);
    m_clockwise = true;
    m_radius = 2.0;
    m_currVel = btVector3(0,0,0);
    m_currAcc = btVector3(0,0,0);
    positions.clear();
    use_constraint = false;
    m_avgAcc = 0.0f;
    m_state = GlobalState::GetInstance();
	//joint_lengths[0] = 0.0828054 * LENGTH_RATIO;
	//joint_lengths[1] = 0.0828054 * LENGTH_RATIO;
	//joint_lengths[2] = 0.0728054 * LENGTH_RATIO;
	//joint_lengths[3] = 0.0828054 * LENGTH_RATIO;
	//joint_lengths[4] = 0.0728054 * LENGTH_RATIO;
	//joint_lengths[5] = 0.0728054 * LENGTH_RATIO;

    joint_lengths[0] = 0.1;
	joint_lengths[1] = 0.1;
	joint_lengths[2] = 0.1;
	joint_lengths[3] = 0.1;
	joint_lengths[4] = 0.1;
	joint_lengths[5] = 0.1;
	linkHalfExtents = btVector3(1, 0.05 / 2.0 * LENGTH_RATIO, 0.01 / 2.0 * LENGTH_RATIO);
}

Skeleton::~Skeleton()
{
}

void Skeleton::initPhysics()
{
    int upAxis = 1;

    m_guiHelper->setUpAxis(upAxis);

    // set btDefaultCollisionConfiguration, dispatcher, solver and other things
    this->createEmptyDynamicsWorld(m_solverType);

    bool isPreTick = false;
    m_dynamicsWorld->setInternalTickCallback(OnInternalTickCallback, this, isPreTick);

    m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);
    if (m_dynamicsWorld->getDebugDrawer())
    {
        int mode = btIDebugDraw::DBG_DrawWireframe + btIDebugDraw::DBG_DrawContactPoints + btIDebugDraw::DBG_DrawAabb;
        m_dynamicsWorld->getDebugDrawer()->setDebugMode(mode);  //+btIDebugDraw::DBG_DrawConstraintLimits);
    }

    const bool floating = true;
    const bool damping = true;   // disable bullet internal damp
    const bool gyro = false;
    const bool canSleep = false;
    const bool selfCollide = false;

    //std::tuple<btVector3, btScalar> r = m_move.getPos(0, 0.041666667);
	//btVector3 init_pos = btVector3(0.71554315598038409, 0.0000000000000000, -0.23249407030114064) * LENGTH_RATIO;
	//btScalar omega = std::get<1>(r);
    //btVector3 init_pos(0,0,-0.0752367);
//    btVector3 init_pos(0,0,-0.08);
	btVector3 init_pos(0, 0, 0);

    /////////////////////////////////////////////////////////////////
    // construct the main skeleton
    /////////////////////////////////////////////////////////////////
	{
        btVector3 baseInertiaDiag(0.f, 0.f, 0.f);
		float baseMass = 0.01f * MASS_RATIO;

		if (baseMass)
		{
			btCollisionShape* pTempBox = new btSphereShape(btScalar(0.01));
			pTempBox->calculateLocalInertia(baseMass, baseInertiaDiag);
			delete pTempBox;
		}

		customMultiBody* pMultiBody = new customMultiBody(m_numLinks, baseMass, baseInertiaDiag, !floating, canSleep);
		//pMultiBody->useRK4Integration(true);

		// set base position
		btScalar angle;
		angle = 0 * SIMD_PI / 180.f;
		btQuaternion baseQ = btQuaternion(btVector3(0, 1, 0).normalized(), angle);
		//    btQuaternion baseQ(0.f, 0.f, 0.f, 1.f);
		pMultiBody->setWorldToBaseRot(baseQ);
		//    btVector3 basePos = btVector3(0.0, 0.0, m_radius);
		pMultiBody->setBasePos(init_pos);

		m_multiBody = pMultiBody;

		for (int i = 0; i < m_numLinks; ++i)
		{
			float linkMass = 1.f * MASS_RATIO;
			btVector3 linkInertiaDiag(0.f, 0.f, 0.f);
			linkHalfExtents[0] = joint_lengths[i] / 2;
			btCollisionShape* shape = 0;
			{
				shape = new btBoxShape(linkHalfExtents);
			}
			shape->calculateLocalInertia(linkMass, linkInertiaDiag);

			// change the linkInertiaDiag to be like a sphere
			float inertia_max = 0.0f;
			if (linkInertiaDiag.x() > inertia_max)
				inertia_max = linkInertiaDiag.x();
			if (linkInertiaDiag.y() > inertia_max)
				inertia_max = linkInertiaDiag.y();
			if (linkInertiaDiag.z() > inertia_max)
				inertia_max = linkInertiaDiag.z();
			linkInertiaDiag = btVector3(inertia_max, inertia_max, inertia_max);
			printf("joint: %d linkInertiaDiag: %f, %f, %f\n", i, linkInertiaDiag.x(), linkInertiaDiag.y(), linkInertiaDiag.z());

			delete shape;

			btVector3 temp;
			btScalar prev_half_length = (i == 0) ? 0 : joint_lengths[i - 1] / 2;
			temp = btVector3(1, 0, 0) * prev_half_length;
			btVector3 parentComToCurrentCom(temp);

			btScalar curr_length = joint_lengths[i];
			temp = btVector3(1, 0, 0) * curr_length / 2;
			btVector3 currentPivotToCurrentCom(temp);

			pMultiBody->setupSpherical(i, 1.0, linkInertiaDiag, i - 1,
									   btQuaternion(0.f, 0.f, 0.f, 1.f),
									   parentComToCurrentCom,
									   currentPivotToCurrentCom, true);
		}

		// init params
		pMultiBody->finalizeMultiDof();
		m_dynamicsWorld->addMultiBody(pMultiBody);
		pMultiBody->setCanSleep(canSleep);
		pMultiBody->setHasSelfCollision(selfCollide);
		pMultiBody->setUseGyroTerm(gyro);
		if (damping)
		{
			btScalar linearDamp = 1.5f;
			btScalar angularDamp = 3.5f;

			// TODO set linear and angular damp for each joint
			pMultiBody->setLinearDamping(linearDamp);
			pMultiBody->setAngularDamping(angularDamp);

			btAlignedObjectArray<btScalar> damps;
			damps.resize(m_numLinks);

			for (int i = 0; i < damps.size(); i++)
			{
				damps[i] = linearDamp;
			}
			pMultiBody->setLinearDampingK1(damps);
			pMultiBody->setLinearDampingK2(damps);

			for (int i = 0; i < damps.size(); i++)
			{
				damps[i] = angularDamp;
			}
			pMultiBody->setAngularDampingK1(damps);
			pMultiBody->setAngularDampingK2(damps);
		}
         
        // init pose
		//{
		//	for (int i = 0; i < pMultiBody->getNumLinks(); i++)
		//	{
		//        if (i == -1) {
		//            btQuaternion q(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
		//			pMultiBody->setJointPosMultiDof(i, q);
		//        } else {
		//			btQuaternion q(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
		//			pMultiBody->setJointPosMultiDof(i, q);
		//        }
		//    }
		//}
		pMultiBody->setJointPosMultiDof(0, btQuaternion(0.270598, 0.270598, -0.653282, 0.653282));
		//pMultiBody->setJointPosMultiDof(0, btQuaternion(0, 0, 0, 1));
		pMultiBody->setJointPosMultiDof(1, btQuaternion(0, 0, 0, 1));
		pMultiBody->setJointPosMultiDof(2, btQuaternion(0, 0, 0, 1));
		pMultiBody->setJointPosMultiDof(3, btQuaternion(0, 0, 0, 1));
		pMultiBody->setJointPosMultiDof(4, btQuaternion(0, 0, 0, 1));
		pMultiBody->setJointPosMultiDof(5, btQuaternion(0, 0, 0, 1));

		// init balance
		m_balanceRot.resize(pMultiBody->getNumLinks());
		m_balanceRot[0] = btQuaternion(0.5, 0.5, -0.5, 0.5);
		//m_balanceRot[0] = btQuaternion(0, 0, 0, 1);
		m_balanceRot[1] = btQuaternion(0, 0, 0, 1);
		m_balanceRot[2] = btQuaternion(0, 0, 0, 1);
		m_balanceRot[3] = btQuaternion(0, 0, 0, 1);
		m_balanceRot[4] = btQuaternion(0, 0, 0, 1);
		m_balanceRot[5] = btQuaternion(0, 0, 0, 1);

		m_Ks.resize(pMultiBody->getNumLinks());
		for (int i = 0; i < m_numLinks; i++)
		{
			m_Ks[i] = 3.5f;
		}

		addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);
	}

    // init p2p constraint
	{
		btVector3 pointInA = m_multiBody->worldPosToLocal(0, init_pos);
		btVector3 pointInB = init_pos;
		m_p2p = new btMultiBodyPoint2Point(m_multiBody, 0, 0, pointInA, pointInB);
		m_p2p->setMaxAppliedImpulse(100);
		m_p2p->setErp(0.8);
		m_dynamicsWorld->addMultiBodyConstraint(m_p2p);
	}

    // change init pos to the tail of the main skeleton
	init_pos = btVector3(0.1 * m_numLinks, 0, 0);

    /////////////////////////////////////////////////////////////////
	// construct the branch skeleton 1
	/////////////////////////////////////////////////////////////////
	//{
	//	int numLink = m_numLinks1;
	//	btVector3 baseInertiaDiag(0.f, 0.f, 0.f);
	//	float baseMass = 0.01f * MASS_RATIO;
	//	if (baseMass)
	//	{
	//		btCollisionShape* pTempBox = new btSphereShape(btScalar(0.01));
	//		pTempBox->calculateLocalInertia(baseMass, baseInertiaDiag);
	//		delete pTempBox;
	//	}
	//	customMultiBody* pMultiBody = new customMultiBody(numLink, baseMass, baseInertiaDiag, !floating, canSleep);
	//	btScalar angle;
	//	angle = 30 * SIMD_PI / 180.f;
	//	btQuaternion baseQ = btQuaternion(btVector3(0, 1, 0).normalized(), angle);
	//	pMultiBody->setWorldToBaseRot(baseQ);
	//	pMultiBody->setBasePos(init_pos);

	//	m_multiBody1 = pMultiBody;

	//	for (int i = 0; i < numLink; ++i)
	//	{
	//		float linkMass = 1.f * MASS_RATIO;
	//		btVector3 linkInertiaDiag(0.f, 0.f, 0.f);
	//		linkHalfExtents[0] = joint_lengths[i] / 2;
	//		btCollisionShape* shape = 0;
	//		{
	//			shape = new btBoxShape(linkHalfExtents);
	//		}
	//		shape->calculateLocalInertia(linkMass, linkInertiaDiag);

	//		// change the linkInertiaDiag to be like a sphere
	//		float inertia_max = 0.0f;
	//		if (linkInertiaDiag.x() > inertia_max)
	//			inertia_max = linkInertiaDiag.x();
	//		if (linkInertiaDiag.y() > inertia_max)
	//			inertia_max = linkInertiaDiag.y();
	//		if (linkInertiaDiag.z() > inertia_max)
	//			inertia_max = linkInertiaDiag.z();
	//		linkInertiaDiag = btVector3(inertia_max, inertia_max, inertia_max);

	//		delete shape;

	//		btVector3 temp;
	//		btScalar prev_half_length = (i == 0) ? 0 : joint_lengths[i - 1] / 2;
	//		temp = btVector3(1, 0, 0) * prev_half_length;
	//		btVector3 parentComToCurrentCom(temp);

	//		btScalar curr_length = joint_lengths[i];
	//		temp = btVector3(1, 0, 0) * curr_length / 2;
	//		btVector3 currentPivotToCurrentCom(temp);

	//		pMultiBody->setupSpherical(i, 1.0, linkInertiaDiag, i - 1,
	//								   btQuaternion(0.f, 0.f, 0.f, 1.f),
	//								   parentComToCurrentCom,
	//								   currentPivotToCurrentCom, true);
	//	}

	//	// init params
	//	pMultiBody->finalizeMultiDof();
	//	m_dynamicsWorld->addMultiBody(pMultiBody);
	//	pMultiBody->setCanSleep(canSleep);
	//	pMultiBody->setHasSelfCollision(selfCollide);
	//	pMultiBody->setUseGyroTerm(gyro);
	//	if (damping)
	//	{
	//		btScalar linearDamp = 1.5f;
	//		btScalar angularDamp = 3.5f;

	//		// TODO set linear and angular damp for each joint
	//		pMultiBody->setLinearDamping(linearDamp);
	//		pMultiBody->setAngularDamping(angularDamp);

	//		btAlignedObjectArray<btScalar> damps;
	//		damps.resize(m_numLinks);

	//		for (int i = 0; i < damps.size(); i++)
	//		{
	//			damps[i] = linearDamp;
	//		}
	//		pMultiBody->setLinearDampingK1(damps);
	//		pMultiBody->setLinearDampingK2(damps);

	//		for (int i = 0; i < damps.size(); i++)
	//		{
	//			damps[i] = angularDamp;
	//		}
	//		pMultiBody->setAngularDampingK1(damps);
	//		pMultiBody->setAngularDampingK2(damps);
	//	}

 //       // init pose
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		pMultiBody->setJointPosMultiDof(i, btQuaternion(0, 0, 0, 1));
	//	}

	//	// init balance
	//	m_balanceRot1.resize(pMultiBody->getNumLinks());
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		m_balanceRot1[i] = btQuaternion(0, 0, 0, 1);
	//	}

 //       // use m_ks from main skeleton

	//	addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);
	//}

 //       // init p2p constraint
	//{
	//	btVector3 pointInA = m_multiBody->worldPosToLocal(m_numLinks - 1, init_pos);
	//	btVector3 pointInB = m_multiBody1->worldPosToLocal(0, init_pos);
	//	btMultiBodyPoint2Point* p2p = new btMultiBodyPoint2Point(m_multiBody, m_numLinks - 1, m_multiBody1, 0, pointInA, pointInB);
	//	p2p->setMaxAppliedImpulse(100);
	//	p2p->setErp(0.8);
	//	m_dynamicsWorld->addMultiBodyConstraint(p2p);
	//}


    /////////////////////////////////////////////////////////////////
	// construct the branch skeleton 2
	/////////////////////////////////////////////////////////////////
	//{
	//	int numLink = m_numLinks2;
	//	btVector3 baseInertiaDiag(0.f, 0.f, 0.f);
	//	float baseMass = 0.01f * MASS_RATIO;
	//	if (baseMass)
	//	{
	//		btCollisionShape* pTempBox = new btSphereShape(btScalar(0.01));
	//		pTempBox->calculateLocalInertia(baseMass, baseInertiaDiag);
	//		delete pTempBox;
	//	}
	//	customMultiBody* pMultiBody = new customMultiBody(numLink, baseMass, baseInertiaDiag, !floating, canSleep);
	//	btScalar angle;
	//	angle = -30 * SIMD_PI / 180.f;
	//	btQuaternion baseQ = btQuaternion(btVector3(0, 1, 0).normalized(), angle);
	//	pMultiBody->setWorldToBaseRot(baseQ);
	//	pMultiBody->setBasePos(init_pos);

	//	m_multiBody2 = pMultiBody;

	//	for (int i = 0; i < numLink; ++i)
	//	{
	//		float linkMass = 1.f * MASS_RATIO;
	//		btVector3 linkInertiaDiag(0.f, 0.f, 0.f);
	//		linkHalfExtents[0] = joint_lengths[i] / 2;
	//		btCollisionShape* shape = 0;
	//		{
	//			shape = new btBoxShape(linkHalfExtents);
	//		}
	//		shape->calculateLocalInertia(linkMass, linkInertiaDiag);

	//		// change the linkInertiaDiag to be like a sphere
	//		float inertia_max = 0.0f;
	//		if (linkInertiaDiag.x() > inertia_max)
	//			inertia_max = linkInertiaDiag.x();
	//		if (linkInertiaDiag.y() > inertia_max)
	//			inertia_max = linkInertiaDiag.y();
	//		if (linkInertiaDiag.z() > inertia_max)
	//			inertia_max = linkInertiaDiag.z();
	//		linkInertiaDiag = btVector3(inertia_max, inertia_max, inertia_max);

	//		delete shape;

	//		btVector3 temp;
	//		btScalar prev_half_length = (i == 0) ? 0 : joint_lengths[i - 1] / 2;
	//		temp = btVector3(1, 0, 0) * prev_half_length;
	//		btVector3 parentComToCurrentCom(temp);

	//		btScalar curr_length = joint_lengths[i];
	//		temp = btVector3(1, 0, 0) * curr_length / 2;
	//		btVector3 currentPivotToCurrentCom(temp);

	//		pMultiBody->setupSpherical(i, 1.0, linkInertiaDiag, i - 1,
	//								   btQuaternion(0.f, 0.f, 0.f, 1.f),
	//								   parentComToCurrentCom,
	//								   currentPivotToCurrentCom, true);
	//	}

	//	// init params
	//	pMultiBody->finalizeMultiDof();
	//	m_dynamicsWorld->addMultiBody(pMultiBody);
	//	pMultiBody->setCanSleep(canSleep);
	//	pMultiBody->setHasSelfCollision(selfCollide);
	//	pMultiBody->setUseGyroTerm(gyro);
	//	if (damping)
	//	{
	//		btScalar linearDamp = 1.5f;
	//		btScalar angularDamp = 3.5f;

	//		// TODO set linear and angular damp for each joint
	//		pMultiBody->setLinearDamping(linearDamp);
	//		pMultiBody->setAngularDamping(angularDamp);

	//		btAlignedObjectArray<btScalar> damps;
	//		damps.resize(m_numLinks);

	//		for (int i = 0; i < damps.size(); i++)
	//		{
	//			damps[i] = linearDamp;
	//		}
	//		pMultiBody->setLinearDampingK1(damps);
	//		pMultiBody->setLinearDampingK2(damps);

	//		for (int i = 0; i < damps.size(); i++)
	//		{
	//			damps[i] = angularDamp;
	//		}
	//		pMultiBody->setAngularDampingK1(damps);
	//		pMultiBody->setAngularDampingK2(damps);
	//	}

	//	// init pose
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		pMultiBody->setJointPosMultiDof(i, btQuaternion(0, 0, 0, 1));
	//	}

	//	// init balance
	//	m_balanceRot2.resize(pMultiBody->getNumLinks());
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		m_balanceRot2[i] = btQuaternion(0, 0, 0, 1);
	//	}

	//	// use m_ks from main skeleton

	//	addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);
	//}

	//// init p2p constraint
	//{
	//	btVector3 pointInA = m_multiBody->worldPosToLocal(m_numLinks - 1, init_pos);
	//	btVector3 pointInB = m_multiBody2->worldPosToLocal(0, init_pos);
	//	btMultiBodyPoint2Point* p2p = new btMultiBodyPoint2Point(m_multiBody, m_numLinks - 1, m_multiBody2, 0, pointInA, pointInB);
	//	p2p->setMaxAppliedImpulse(100);
	//	p2p->setErp(0.8);
	//	m_dynamicsWorld->addMultiBodyConstraint(p2p);
	//}

    /////////////////////////////////////////////////////////////////
    // construct the marker
    /////////////////////////////////////////////////////////////////
    {
		btCollisionShape* colShape = new btSphereShape(btScalar(0.02 * LENGTH_RATIO));

        /// Create Dynamic Objects
        btTransform startTransform;
        startTransform.setIdentity();
        startTransform.setOrigin(init_pos);

        btScalar mass(0.f);

        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);

        btVector3 localInertia(0, 0, 0);
        if (isDynamic)
            colShape->calculateLocalInertia(mass, localInertia);

        //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
        btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
        m_marker = new btRigidBody(rbInfo);
        m_marker->setRestitution(0.0);
        m_dynamicsWorld->addRigidBody(m_marker, NOTHING, NOTHING);

        if (!WIRE_FRAME)
        {
            m_guiHelper->createCollisionShapeGraphicsObject(colShape);
            btVector4 color(1, 0, 0, 1);
            m_guiHelper->createCollisionObjectGraphicsObject(dynamic_cast<btCollisionObject*>(m_marker), color);
        }
    }

    ///////////////////////////////////////////////////////////////////
    //// construct the marker
    ///////////////////////////////////////////////////////////////////
    //{
    //    btCollisionShape* colShape = new btSphereShape(btScalar(0.02));

    //    /// Create Dynamic Objects
    //    btTransform startTransform;
    //    startTransform.setIdentity();
    //    startTransform.setOrigin(btVector3(0, -0.33, -0.405));

    //    btScalar mass(0.f);

    //    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    //    bool isDynamic = (mass != 0.f);

    //    btVector3 localInertia(0, 0, 0);
    //    if (isDynamic)
    //        colShape->calculateLocalInertia(mass, localInertia);

    //    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    //    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    //    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    //    m_marker = new btRigidBody(rbInfo);
    //    m_marker->setRestitution(0.0);
    //    m_dynamicsWorld->addRigidBody(m_marker, NOTHING, NOTHING);

    //    if (!WIRE_FRAME)
    //    {
    //        m_guiHelper->createCollisionShapeGraphicsObject(colShape);
    //        btVector4 color(1, 0, 0, 1);
    //        m_guiHelper->createCollisionObjectGraphicsObject(dynamic_cast<btCollisionObject*>(m_marker), color);
    //    }
    //}

    /////////////////////////////////////////////////////////////////
    // construct the collider
    /////////////////////////////////////////////////////////////////
  //  {
  //      btCollisionShape* colShape = new btCapsuleShape(btScalar(0.06), btScalar(0.18));

  //      /// Create Dynamic Objects
  //      btTransform startTransform;
  //      startTransform.setIdentity();
		//startTransform.setOrigin(btVector3(0.0, -0.2, 0.7523666 * LENGTH_RATIO));

  //      btScalar mass(0.f);

  //      //rigidbody is dynamic if and only if mass is non zero, otherwise static
  //      bool isDynamic = (mass != 0.f);

  //      btVector3 localInertia(0, 0, 0);
  //      if (isDynamic)
  //          colShape->calculateLocalInertia(mass, localInertia);

  //      //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
  //      btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
  //      btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
  //      m_collider = new btRigidBody(rbInfo);
  //      m_collider->setRestitution(0.0);
  //      m_dynamicsWorld->addRigidBody(m_collider, TARGET_BODY, BONE_BODY);

  //      if (!WIRE_FRAME)
  //      {
  //          m_guiHelper->createCollisionShapeGraphicsObject(colShape);
  //          btVector4 color(1, 0, 1, 1);
  //          m_guiHelper->createCollisionObjectGraphicsObject(dynamic_cast<btCollisionObject*>(m_collider), color);
  //      }
  //  }


    m_dynamicsWorld->setGravity(btVector3(0, 0, 0));

    btContactSolverInfo& si = m_dynamicsWorld->getSolverInfo();
    si.m_numIterations = 10;
    si.m_globalCfm = 0.05f;
    si.m_erp2 = 0.1f;

    // setMaxCoordinateVelocity will also limit the base vel and rot. so we use setMaxOmega
//    m_multiBody->setMaxCoordinateVelocity(2.0);
    m_multiBody->setMaxOmega(100.0f);
    // setMaxOmegaX will cause instability problem.
//    m_multiBody->setMaxOmegaX(0.5);
}

std::pair<btVector3, btVector3> Skeleton::applySpringForce(customMultiBody* mb, const btAlignedObjectArray<btQuaternion>& balanceRot, const btVector3& tailTorque)
{
    // if reduce_factor > 0.8, it looks unnatural. 
	btScalar reduce_factor = 0.6;

    btAlignedObjectArray<btVector3> fbtorques;
	fbtorques.resize(mb->getNumLinks(), btVector3(0,0,0));

    const btScalar max_torque = 0.1f * MASS_RATIO;
    const btVector3 joint_dir = btVector3(1, 0, 0);
	btQuaternion prevQ = mb->getWorldToBaseRot().inverse();

    const btVector3 joint_y_dir = btVector3(0, 1, 0);
	btQuaternion prevXQ = mb->getWorldToBaseRot().inverse();

    // forwards
    for (int i = 0; i < mb->getNumLinks(); ++i)
    {
		btVector3 c = mb->getLink(i).m_cachedWorldTransform.getOrigin();
        // cause joint_dir is a unit vector, so after localDirToWorld it does not need normalization.
		btVector3 dir = mb->localDirToWorld(i, joint_dir);
        btVector3 origin = c + (-dir) * joint_lengths[i] * 0.5;
		btQuaternion balanceQ = balanceRot[i];

        {
			// reverse x axis rot
			btQuaternion currentQuat(mb->getJointPosMultiDof(i)[0],
									 mb->getJointPosMultiDof(i)[1],
									 mb->getJointPosMultiDof(i)[2],
									 mb->getJointPosMultiDof(i)[3]);

            btVector3 angleDiff;
            btQuaternion cq = prevXQ * currentQuat;
			btQuaternion tq = prevXQ * balanceQ;
			btQuaternion middle = cq.slerp(tq, 0.5);
			btQuaternion relRot = middle.inverse() * tq;
			btGeneric6DofSpring2Constraint::matrixToEulerXYZ(btMatrix3x3(relRot), angleDiff);
			printf("joint %d : middle diff: %f, %f, %f\n", i, rad2degree(angleDiff.x()), rad2degree(angleDiff.y()), rad2degree(angleDiff.z()));

            btVector3 torque;
            if (fabs(angleDiff.x()) > 1 * SIMD_PI / 180.f)
			{
				btMatrix3x3 Rt = btMatrix3x3(cq);
				btMatrix3x3 Ibody;
				btVector3 inertia = mb->getLink(i).m_inertiaLocal;
				Ibody.setValue(inertia.x(), 0, 0,
							   0, inertia.y(), 0,
							   0, 0, inertia.z());
				btMatrix3x3 Iinv = Rt * Ibody.inverse() * Rt.transpose();
				btVector3 Ivx = Iinv * btVector3(1, 0, 0);
				btScalar Ix = Ivx.x();
				btScalar kw = 0.006;
				btScalar phi = angleDiff.x() / SIMD_PI * 2.0;
				btScalar t = exp(-phi * phi);
				btScalar scaling = 5.0;
				btScalar ratio = scaling * kw * exp(-btSqrt(kw * Ix) * t);
				btScalar positive = btIsNegative(angleDiff.x()) ? -1 : 1;
				torque = -dir * positive * ratio;
                    
				//btScalar ratio = pow(fabs(angleDiff.x() / SIMD_PI), 2);
				//torque = -dir * angleDiff.x() * ratio * 5;
				if (torque.safeNorm() > max_torque)
				{
					torque *= max_torque / torque.safeNorm();
				}
				mb->addLinkTorque(i, torque);
				
                fbtorques[i] += -torque;
			}

            if (WIRE_FRAME)
			{
				m_guiHelper->getRenderInterface()->drawLine(c, c + torque.safeNormalize() * 0.05, btVector4(0, 0, 1, 0), btScalar(2));
			}
			prevXQ = prevXQ * currentQuat;
		}

        {
			btVector3 temp1 = quatRotate(prevQ * balanceQ, joint_dir * joint_lengths[i]);
			btVector3 balanceP = temp1 * 0.5 + origin;

			temp1 = balanceP - origin;
			btVector3 temp2 = temp1.normalized();
			float angle = btAcos(temp2.dot(dir));
			btVector3 dist2balance = balanceP - c - (balanceP - c).dot(dir) * dir;

			if (angle > 1 * SIMD_PI / 180.f)
			{
				btVector3 torque = dir.cross(dist2balance.normalized()) * m_Ks[i] * joint_lengths[i] * angle;
				if (torque.safeNorm() > max_torque)
				{
					torque *= max_torque / torque.safeNorm();
				}
				mb->addLinkTorque(i, torque);
				//printf("joint: %d : torque %f, %f, %f\n", i, torque.x(), torque.y(), torque.z());

				fbtorques[i] += -torque;
			}
			if (WIRE_FRAME)
			{
				//m_guiHelper->getRenderInterface()->drawLine(c, c + torque, btVector4(0, 0, 1, 0), btScalar(2));
			}
			prevQ = prevQ * balanceQ;
		}
    }

    // backwards
    for (int i = mb->getNumLinks() - 1; i >= 0; --i)
	{
		if (i == mb->getNumLinks() - 1)
			mb->addLinkTorque(i, tailTorque * reduce_factor);
		else
		{
			mb->addLinkTorque(i, fbtorques[i + 1] * reduce_factor);
        }
	}
    
    return std::make_pair(fbtorques[0], -fbtorques[mb->getNumLinks() - 1]);
}

void Skeleton::addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const float* joint_lengths)
{
    btAlignedObjectArray<btQuaternion> world_to_local;
    world_to_local.resize(pMultiBody->getNumLinks() + 1);

    btAlignedObjectArray<btVector3> local_origin;
    local_origin.resize(pMultiBody->getNumLinks() + 1);
    world_to_local[0] = pMultiBody->getWorldToBaseRot();
    local_origin[0] = pMultiBody->getBasePos();

    for (int i = 0; i < pMultiBody->getNumLinks(); ++i)
    {
        const int parent = pMultiBody->getParent(i);
        world_to_local[i + 1] = pMultiBody->getParentToLocalRot(i) * world_to_local[parent + 1];
        local_origin[i + 1] = local_origin[parent + 1] + (quatRotate(world_to_local[i + 1].inverse(), pMultiBody->getRVector(i)));
    }

    for (int i = 0; i < pMultiBody->getNumLinks(); ++i)
    {
        btVector3 posr = local_origin[i + 1];
        btScalar quat[4] = {-world_to_local[i + 1].x(), -world_to_local[i + 1].y(), -world_to_local[i + 1].z(), world_to_local[i + 1].w()};

        linkHalfExtents[0] = joint_lengths[i] / 2.0f;
		btCollisionShape* shape = new btBoxShape(linkHalfExtents);
        btMultiBodyLinkCollider* col = new btMultiBodyLinkCollider(pMultiBody, i);
        col->setCollisionShape(shape);

        btTransform tr;
        tr.setIdentity();
        tr.setOrigin(posr);
        tr.setRotation(btQuaternion(quat[0], quat[1], quat[2], quat[3]));
        col->setWorldTransform(tr);

        col->setRestitution(0.0);

        // set collision group and mask, only collide with objects with mask is true where it collide with collider
        m_dynamicsWorld->addCollisionObject(col, BONE_BODY, TARGET_BODY);  //,2,1+2);

        if (!WIRE_FRAME)
        {
            m_guiHelper->createCollisionShapeGraphicsObject(shape);
            btVector4 color(0, 1, 0, 1);
            m_guiHelper->createCollisionObjectGraphicsObject(col, color);
        }

        pMultiBody->getLink(i).m_collider = col;
    }
}

void Skeleton::moveCollider(const btVector3& pos){
    if ( m_collider )
    {
        btTransform tr;
        tr.setIdentity();
        tr.setOrigin(pos);
        tr.setRotation(btQuaternion(0.0, 0.0, 0.0, 1.0));
        m_collider->setWorldTransform(tr);
    }
}

void Skeleton::moveMarker(const btVector3& pos){
    if ( m_marker )
    {
        btTransform tr;
        tr.setIdentity();
        tr.setOrigin(pos);
        tr.setRotation(btQuaternion(0.0, 0.0, 0.0, 1.0));
        m_marker->setWorldTransform(tr);
    }
}

void Skeleton::getLinearAcc(float deltaTime, int m_step, const btVector3& pos)
{
    if ( positions.size() == POS_ARRAY_SIZE)
    {
        positions.pop_front();
    }
    positions.push_back(pos);

    if ( positions.size() > 1 ) {
        int length = positions.size();
        m_currVel = (positions[length - 1] - positions[length - 2]) / deltaTime;
    }

    if ( positions.size() > 2 ) {
        int length = positions.size();
        m_currAcc = (positions[length - 1] + positions[length - 3] - 2 * positions[length - 2]) / deltaTime / deltaTime;

        if ( accs.size() == ACC_ARRAY_SIZE)
        {
            accs.pop_front();
        }
        accs.push_back(m_currAcc);
    }
//    printf("acc %f\n", m_currAcc.norm());
//    printf("vel %f\n", m_currVel.norm());
}

void Skeleton::applyBaseLinearDragForce(const btVector3& dir)
{
    btVector3 currBaseVel;
    float scaling = 0.4;
    float damp = 0.0; // if damp is set to 0.2 or more, the tail will wind up in horizontal oscillating motion

    // TODO replace to avg acc
    if ( m_avgAcc > 0.1 ) {
        btVector3 g = dir * m_avgAcc * scaling;
        for (int i = 0; i < m_numLinks; ++i) {
            m_multiBody->addLinkForce(i, g * m_multiBody->getLink(i).m_mass / pow(i+1, damp));
        }
    }
}

void Skeleton::applyBaseCentrifugalForce(float deltaTime, const btVector3& pos, float omega)
{
    btVector3 dir = pos - btVector3(0,0,0);

    std::set<int> links = m_state->getCollisionEvent();
    int min_link = m_multiBody->getNumLinks();
    for ( auto iter = links.begin(); iter != links.end(); iter++ )
    {
        if (min_link > *iter)
            min_link = *iter;
    }
    min_link += 1;
    if (min_link >= m_multiBody->getNumLinks())
        min_link = m_multiBody->getNumLinks();

    const btScalar scaling = 1.0;
    for ( int i = 0; i < min_link; i++ )
    {
        btVector3 c = m_multiBody->getLink(i).m_cachedWorldTransform.getOrigin();
        btScalar r = sqrt(c.x() * c.x() + c.z() * c.z());
        btVector3 force = m_multiBody->getLinkMass(i) * omega * omega * r * dir.normalized() * (1 / pow(i + 1, 0.6));
        m_multiBody->addLinkForce(i, force);

        if ( WIRE_FRAME ) {
            m_guiHelper->getRenderInterface()->drawLine(c, c + force, btVector4(1, 1, 0, 0),
                                                        btScalar(2));
        }
    }
}

void Skeleton::limitMaxTwist(float max_angle)
{
    for (int i = 1; i < m_multiBody->getNumLinks(); ++i) {
        btQuaternion q(m_multiBody->getJointPosMultiDof(i)[0],
                       m_multiBody->getJointPosMultiDof(i)[1],
                       m_multiBody->getJointPosMultiDof(i)[2],
                       m_multiBody->getJointPosMultiDof(i)[3]);

        bool needed = false;
        float a = q.getAngle();
        if ( a > max_angle ) {
            a = max_angle;
            needed = true;
        }
        if ( a < -max_angle) {
            a = -max_angle;
            needed = true;
        }
        if (needed) {

            btQuaternion newq;
            if ( fabs(sin(q.getAngle() / 2 )) > 1e-4 )
            {
                newq[3] = cos(a / 2);
                float ratio = sin(q.getAngle() / 2) / sin(q.getAngle() / 2);
                newq[0] = q[0] * ratio;
                newq[1] = q[1] * ratio;
                newq[2] = q[2] * ratio;
            }
            m_multiBody->setJointPosMultiDof(i, &newq[0]);
        }
    }
}

void Skeleton::OnInternalTickCallback(btDynamicsWorld* world, btScalar timeStep)
{
    Skeleton* demo = static_cast<Skeleton*>(world->getWorldUserInfo());

    btVector3 feedbackTorque;
	//feedbackTorque = demo->applySpringForce(demo->m_multiBody1, demo->m_balanceRot1);
    demo->applySpringForce(demo->m_multiBody, demo->m_balanceRot);
	
	//demo->applySpringForce(demo->m_multiBody2, demo->m_balanceRot2);

    // draw center of the link
    if (WIRE_FRAME)
    {
        btVector4 color(1, 0, 0, 1);
        btVector3 base = demo->m_multiBody->getBasePos();
        for ( int i = 0; i < demo->m_multiBody->getNumLinks(); i++ ){
            btVector3 c = demo->m_multiBody->getLink(i).m_cachedWorldTransform.getOrigin();
            demo->m_guiHelper->getRenderInterface()->drawPoint(c, color, btScalar(5));
        }
    }

    int numManifolds = world->getDispatcher()->getNumManifolds();
    bool has_collision = false;
    for (int i = 0; i < numManifolds; i++)
    {
        btPersistentManifold* contactManifold = world->getDispatcher()->getManifoldByIndexInternal(i);
        const btCollisionObject* obA = contactManifold->getBody0();
        const btCollisionObject* obB = contactManifold->getBody1();

        // only detect collision between bones and collider
        if (
                (!(obA->getBroadphaseHandle()->m_collisionFilterGroup & BONE_BODY)
                    && !(obB->getBroadphaseHandle()->m_collisionFilterGroup | BONE_BODY))
                ||
                (!(obA->getBroadphaseHandle()->m_collisionFilterGroup & TARGET_BODY)
                 && !(obB->getBroadphaseHandle()->m_collisionFilterGroup | TARGET_BODY))
        )
            continue;

        int numContacts = contactManifold->getNumContacts();
        for (int j = 0; j < numContacts; j++)
        {
            // when there is no collision, set gravity to go downside
            has_collision = true;

            const btCollisionObject* bone;
            if (obA->getBroadphaseHandle()->m_collisionFilterGroup & BONE_BODY)
                bone = obA;
            else
                bone = obB;

            btMultiBodyLinkCollider* link = (btMultiBodyLinkCollider*)btMultiBodyLinkCollider::upcast(bone);
            demo->m_state->insertCollisionEvent(link->m_link);

            btManifoldPoint& pt = contactManifold->getContactPoint(j);
            printf("step : %d, collision impulse: %f\n", demo->m_step, pt.getAppliedImpulse());

            btVector3 ptOnBone, normalOnBone;
            // if obA is bone
            if ( obA->getBroadphaseHandle()->m_collisionFilterGroup ){
                ptOnBone = pt.getPositionWorldOnA();
                normalOnBone = pt.m_normalWorldOnB;
            }
            else{
                ptOnBone = pt.getPositionWorldOnB();
                normalOnBone = -pt.m_normalWorldOnB;
            }
            normalOnBone.normalize();
            printf("ptA: %f, %f, %f\n", ptOnBone.getX(), ptOnBone.getY(), ptOnBone.getZ());
            printf("normal: %f, %f, %f\n", normalOnBone.getX(), normalOnBone.getY(), normalOnBone.getZ());

//            btVector3 to = ptOnBone + normalOnBone * (pt.getAppliedImpulse() + 1e-2) * 10;
            btVector3 to = ptOnBone + normalOnBone;
            btVector4 color(1, 0, 0, 1);
            if (demo->m_guiHelper->getRenderInterface())
            {
                if (WIRE_FRAME)
                {
                    demo->m_guiHelper->getRenderInterface()->drawLine(ptOnBone, to, color, btScalar(1));
                }
            }
        }
    }

    if ( has_collision )
        demo->m_g.startGravity();
    else
        demo->m_g.stopGravity(); // when there is no collision, set gravity to zero
}

void Skeleton::stepSimulation(float deltaTime) {
//    float g = m_g.getGVal();
//    m_dynamicsWorld->setGravity(btVector3(0, g, 0));

    printf("step: %d\n", m_step);

 //   //auto [ basePos, omega ] = m_move.getPos(m_time, deltaTime);
	//std::tuple<btVector3, btScalar> r = m_move.getPos(m_time, deltaTime);
	//btVector3 basePos = std::get<0>(r);
	//btScalar omega = std::get<1>(r);
	////
 //   getLinearAcc(deltaTime, m_step, basePos);

 //   if ( accs.size() == ACC_ARRAY_SIZE ) {
 //       float sum = 0.0f;
 //       for ( int i = 0; i < accs.size(); i++ )
 //           sum += accs[i].norm();
 //       m_avgAcc = sum / ACC_ARRAY_SIZE;
 //   }

 //   btVector3 dir = btVector3(0, -1, 0);

 //   applyBaseLinearDragForce(dir);

//    applyBaseCentrifugalForce(deltaTime, basePos, omega);

    // p2p
    //m_p2p->setPivotInB(basePos);

    //moveMarker(basePos);

//    moveCollider(btVector3(0, -2, 2));

    // capture the frames
	//if (m_step % 10 == 0)
 //   {
        const char* gPngFileName = "multibody";
        sprintf(mfileName, "%s_%d.png", gPngFileName, m_step);
        this->m_guiHelper->getAppInterface()->dumpNextFrameToPng(mfileName);
 //   }

//    m_state->clearCollisionEvent();

    // step and update positions
    int substeps = 10;
    btScalar  fixedTimeStep = deltaTime / substeps;
    m_dynamicsWorld->stepSimulation(deltaTime, substeps, fixedTimeStep);
    if ( WIRE_FRAME ) {
        m_dynamicsWorld->debugDrawWorld();
    }

    // set limit will cause instability when collisions happens
//    btScalar max_angle = 120 * SIMD_PI / 180.f;
//    limitMaxTwist(max_angle);

//    std::set<int> collisionLinks = m_state->getCollisionEvent();
//    for ( auto iter = collisionLinks.begin(); iter != collisionLinks.end(); iter++ )
//    {
//        printf("collision link : %d\n", *iter);
//    }

    m_time += deltaTime;
    m_step += 1;

    if ( m_step == 20)
    {
        printf("step\n");
    }
}


class CommonExampleInterface* Dof6Spring2CreateFunc(CommonExampleOptions& options)
{
//    return new Dof6Spring2Setup(options.m_guiHelper);
	return new Skeleton(options.m_guiHelper);
}
