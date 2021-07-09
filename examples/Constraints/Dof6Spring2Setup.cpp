#include <math.h>
#include <limits>
#include <deque>
#include <tuple>
#include <mutex>
#include <set>
#include <list>
#include <deque>
#include <numeric>

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

// damping functions
class dampingFunc
{
public:
	dampingFunc() {}
	virtual ~dampingFunc() {}
	virtual btScalar getChange(btScalar curr, btScalar target) = 0;
	virtual btScalar getChange(btScalar diff) = 0;
};

class linearDamp : public dampingFunc
{
public:
	linearDamp(btScalar s = 1.0) : m_scaling(s) {}
	virtual ~linearDamp() {}
	virtual btScalar getChange(btScalar curr, btScalar target) { return (target - curr); }
	virtual btScalar getChange(btScalar diff) { return diff * m_scaling; }

private:
	btScalar m_scaling;
};

class powerDamp : public dampingFunc
{
public:
	powerDamp(btScalar s, btScalar p) : m_scaling(s), m_power(p) {}
	virtual ~powerDamp() {}
	virtual btScalar getChange(btScalar curr, btScalar target) { 
		return getChange(target - curr); 
	}
	virtual btScalar getChange(btScalar diff) {
		btClamp(diff, -SIMD_PI, SIMD_PI);
		btScalar ratio = pow(fabs(diff / SIMD_PI), m_power);
		return diff * ratio * m_scaling;
	}

private:
	btScalar m_scaling;
	btScalar m_power;
};

class springDamp : public dampingFunc
{
public:
	springDamp(const btVector3& inertia, btScalar s, btScalar k) : m_scaling(s), m_kw(k)
	{
		m_ibody.setValue(inertia.x(), 0, 0, 0, inertia.y(), 0, 0, 0, inertia.z());
		m_useXAxis = true;
		m_rt.setIdentity();
	}
	virtual ~springDamp() {}
	virtual btScalar getChange(btScalar curr, btScalar target) { return getChange(target - curr); }
	virtual btScalar getChange(btScalar diff) {
		btMatrix3x3 Iinv = m_rt * m_ibody.inverse() * m_rt.transpose();
		btScalar I;
		btScalar phi = diff / SIMD_PI * 2.0;
		btScalar t = exp(-phi * phi);
		btScalar ratio;
		if (m_useXAxis)  // use x axis
		{
			btVector3 Ivx = Iinv * btVector3(1, 0, 0);
			I = Ivx.x();
			ratio = m_kw * exp(-btSqrt(m_kw * I) * t);
			return diff * ratio * m_scaling;
		}
		else
		{
			// use y and z axis, get the most
			btVector3 Ivy = Iinv * btVector3(0, 1, 0);
			I = Ivy.y();
			btScalar ratioy = m_kw * exp(-btSqrt(m_kw * I) * t);

			btVector3 Ivz = Iinv * btVector3(0, 0, 1);
			I = Ivz.z();
			btScalar ratioz = m_kw * exp(-btSqrt(m_kw * I) * t);

			ratio = btMax(ratioy, ratioz);
			return diff * ratio * m_scaling;
		} 
	}
	void setRt(const btQuaternion& cq) {
		m_rt = btMatrix3x3(cq);
	}
	void setXaxis(bool useX){
		m_useXAxis = useX;
	}

private:
	btMatrix3x3 m_ibody;
	btMatrix3x3 m_rt;
	bool m_useXAxis;
	btScalar m_scaling;
	btScalar m_kw;
};

class pdDamp : public dampingFunc
{
public:
	pdDamp(double dt, double max, double Kp, double Kd, double Ki, double e) 
		: _dt(dt),
		  _max(max),
		  _Kp(Kp),
		  _Kd(Kd),
		  _Ki(Ki),
		  _maxError(e),
		  _pre_error(0),
		  _integral(0),
		  _seqsize(5),
		  _isfirst(true)
	{}

	virtual btScalar getChange(btScalar curr, btScalar target) { return getChange(target - curr); }
	virtual btScalar getChange(btScalar error) {
		// Proportional term
		double Pout = _Kp * error;

		// Integral term
		_integral = std::accumulate(_seqs.begin(), _seqs.end(), decltype(_seqs)::value_type(0));
		double Iout = _Ki * (_integral * _dt);
		if (_seqs.size() == _seqsize)
		{
			_seqs.pop_front();
		}
		_seqs.push_back(error);

		// Derivative term
		double Dout;
		if (_isfirst)
		{
			_isfirst = false;
			Dout = 0;
		}
		else
		{
			double derivative = (error - _pre_error) / _dt;
			Dout = _Kd * derivative;
			if (Dout > _maxError)
				Dout = _maxError;
			else if (Dout < -_maxError)
				Dout = _maxError;
		}
		_pre_error = error;

		// Calculate total output
		double output = Pout + Iout + Dout;

		// Restrict to max/min
		if (output > _max)
			output = _max;
		else if (output < -_max)
			output = _max;

		return output;
	}

private:
	double _dt;
	double _max;
	double _min;
	double _Kp;
	double _Kd;
	double _Ki;
	double _pre_error;
	double _integral;
	double _maxError;
	std::deque<double> _seqs;
	int _seqsize;
	bool _isfirst;
};

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
		m_parent = nullptr;
		m_tailPos = btVector3(0, 0, 0);
		m_tailRot = btQuaternion(0, 0, 0, 1);
		m_tailBalanceRot = btQuaternion(0, 0, 0, 1);
		m_fbRadialTorque = btVector3(0, 0, 0);
		m_fbXaxisTorque = btVector3(0, 0, 0); 
		m_firstspringcalc = true;
    }

    // for multi branch
	void setParentCustomBody(customMultiBody * p) { m_parent = p; }
	customMultiBody* getParentCustomBody() { return m_parent; }
	void addChild(customMultiBody* c) { m_childs.push_back(c); }
	std::vector<customMultiBody*> getChilds() { return m_childs; }

    void constructTreeList()
	{
		std::list<customMultiBody*> todos;
		todos.push_back(this);
		while (!todos.empty())
		{
			customMultiBody* curr = todos.front();
			m_tree.push_back(curr);
			for (int i = 0; i < curr->m_childs.size(); i++)
			{
				todos.push_back(curr->m_childs[i]);
            }
			todos.pop_front();
        }
    }

    void updateTailRots()
	{
		for (int i = 0; i < m_tree.size(); i++)
		{
			btQuaternion baseRot = btQuaternion(0, 0, 0, 1);
			btQuaternion baseBalanceRot = btQuaternion(0, 0, 0, 1);
			customMultiBody* p = m_tree[i]->getParentCustomBody();
			if (p != nullptr) {
				baseRot = p->m_tailRot;
				baseBalanceRot = p->m_tailBalanceRot;
			}
			m_tree[i]->updateCurrTailRot(baseRot, baseBalanceRot);
        }
    }

    void updateCurrTailRot(const btQuaternion& prevBaseQuat, const btQuaternion& prevBalanceQuat = btQuaternion(0,0,0,1))
	{
		btQuaternion baseQ = prevBaseQuat;
		for (int i = 0; i < getNumLinks(); i++)
		{
			btQuaternion currentQuat(getJointPosMultiDof(i)[0], getJointPosMultiDof(i)[1], getJointPosMultiDof(i)[2], getJointPosMultiDof(i)[3]);
			baseQ = baseQ * currentQuat;
		}
		m_tailRot = baseQ;

        btQuaternion balanceQ = prevBalanceQuat;
		for (int i = 0; i < getNumLinks(); i++)
		{
			balanceQ = balanceQ * m_balanceRot[i];
		}
		m_tailBalanceRot = balanceQ;
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

public:
	btVector3 m_tailPos;
	btQuaternion m_tailRot;
	btQuaternion m_tailBalanceRot;
	std::deque<customMultiBody*> m_tree;
	btVector3 m_fbRadialTorque;
	btVector3 m_fbXaxisTorque;
	btAlignedObjectArray<btQuaternion> m_balanceRot;
	std::vector<btVector3> m_rtdirs;
	std::vector<int> m_signs;
	bool m_firstspringcalc;

protected:
    btScalar m_maxOmegaX;
    btScalar m_maxOmega;
	customMultiBody* m_parent;
	std::vector<customMultiBody*> m_childs;
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
    float joint_lengths[20];
	btVector3 linkHalfExtents;

public:
    Skeleton(struct GUIHelperInterface* helper);
    virtual ~Skeleton();
    virtual void initPhysics();
    virtual void stepSimulation(float deltaTime);
    virtual void resetCamera()
    {
		float dist = 3.0 * LENGTH_RATIO;
        float pitch = -40;
        float yaw = -135;
		float targetPos[3] = {0.5, 0, 0.0};
        m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
    }
    void addColliders(customMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const float* joint_lengths);
    void moveCollider(const btVector3& pos);
    void moveMarker(const btVector3& pos);
    void applyBaseLinearDragForce(const btVector3& dir);
    void applyBaseCentrifugalForce(float deltaTime, const btVector3& pos, float omega);
    static void OnInternalTickCallback(btDynamicsWorld* world, btScalar timeStep);
    void getLinearAcc(float deltaTime, int m_step, const btVector3& pos);
    void limitMaxTwist(float max_angle);
	void applySpringForce(customMultiBody* mb);
};

Skeleton::Skeleton(struct GUIHelperInterface* helper)
        : CommonMultiBodyBase(helper)
{
    m_time = btScalar(0.0);
    m_step = 1;
    m_numLinks = 5;
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

    for (int i = 0; i < m_numLinks; i++)
	{
		joint_lengths[i] = 0.1;
    }
	linkHalfExtents = btVector3(1, 0.01 / 2.0 * LENGTH_RATIO, 0.05 / 2.0 * LENGTH_RATIO);
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
		btQuaternion baseQ = btQuaternion(0,0,0,1);
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
		for (int i = 0; i < pMultiBody->getNumLinks(); i++)
		{
			if (i == 0)
			{
				btQuaternion q(btVector3(0, 0, 1).normalized(), -86 * SIMD_PI / 180.f);
				pMultiBody->setJointPosMultiDof(i, q);
			}
			else
			{
				btQuaternion q(btVector3(0, 0, 1).normalized(), 0 * SIMD_PI / 180.f);
				pMultiBody->setJointPosMultiDof(i, q);
			}
		}
		// init balance
		pMultiBody->m_balanceRot.resize(pMultiBody->getNumLinks());
		for (int i = 0; i < pMultiBody->getNumLinks(); i++)
		{
			if (i == 0)
			{
				btQuaternion q(btVector3(0, 0, 1).normalized(), -90 * SIMD_PI / 180.f);
				pMultiBody->m_balanceRot[i] = q;
			}
			else
			{
				btQuaternion q(btVector3(0, 0, 1).normalized(), 0 * SIMD_PI / 180.f);
				pMultiBody->m_balanceRot[i] = q;
			}
		}

  //      // init pose
		//for (int i = 0; i < pMultiBody->getNumLinks(); i++)
		//{
		//	btQuaternion q(btVector3(1, 0, 0).normalized(), 0 * SIMD_PI / 180.f);
		//	pMultiBody->setJointPosMultiDof(i, q);
		//}

		//// init balance
		//pMultiBody->m_balanceRot.resize(pMultiBody->getNumLinks());
		//for (int i = 0; i < pMultiBody->getNumLinks(); i++)
		//{
		//	btQuaternion q(btVector3(1, 0, 0).normalized(), 0 * SIMD_PI / 180.f);
		//	pMultiBody->m_balanceRot[i] = q;
		//}

		m_Ks.resize(pMultiBody->getNumLinks());
		for (int i = 0; i < m_numLinks; i++)
		{
			m_Ks[i] = 1.0f;
		}

		addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);

        m_multiBody->updateCurrTailRot(baseQ);
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
	//	m_multiBody1 = pMultiBody;

 //       // construct the multibranch tree
	//	m_multiBody->addChild(m_multiBody1);
	//	m_multiBody1->setParentCustomBody(m_multiBody);

 //       btQuaternion baseQ = btQuaternion(0, 0, 0, 1);
	//	pMultiBody->setWorldToBaseRot(baseQ);
 //       if (m_multiBody1->getParentCustomBody() != nullptr)
	//	{
	//		customMultiBody* parent = m_multiBody1->getParentCustomBody();
	//		pMultiBody->setBasePos(parent->m_tailPos);
 //       }
	//	else
	//	{
	//		pMultiBody->setBasePos(btVector3(0, 0, 0));
 //       }

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
	//		if (i == 0 && m_multiBody1->getParentCustomBody() != nullptr)
	//		{
	//			btQuaternion q = btQuaternion(btVector3(0, 1, 0).normalized(), -30 * SIMD_PI / 180.f);
	//			//btQuaternion q = btQuaternion(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
	//			//btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			btQuaternion tail = m_multiBody1->getParentCustomBody()->m_tailRot;
	//			q = tail * q;
	//			pMultiBody->setJointPosMultiDof(i, q);
 //           }
	//		else
	//		{
	//			btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			pMultiBody->setJointPosMultiDof(i, q);
	//		}
	//	}
	//	//for (int i = 0; i < pMultiBody->getNumLinks(); i++)
	//	//{
	//	//	btQuaternion q(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
	//	//	pMultiBody->setJointPosMultiDof(i, q);
	//	//}

	//	// init balance
	//	pMultiBody->m_balanceRot.resize(pMultiBody->getNumLinks());
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		
	//		if (i == 0 && m_multiBody1->getParentCustomBody() != nullptr)
	//		{
	//			btQuaternion q = btQuaternion(btVector3(0, 1, 0).normalized(), -30 * SIMD_PI / 180.f);
	//			pMultiBody->m_balanceRot[i] = q;
 //           }
	//		else
	//		{
	//			btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			pMultiBody->m_balanceRot[i] = q;
	//		}
	//	}

 //       // use m_ks from main skeleton

	//	addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);
	//}

 //       // init p2p constraint
	//{
	//	btVector3 pointInA = m_multiBody->worldPosToLocal(m_numLinks - 1, m_multiBody->m_tailPos);
	//	btVector3 pointInB = m_multiBody1->worldPosToLocal(0, m_multiBody->m_tailPos);
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
	//	m_multiBody2 = pMultiBody;

	//	// construct the multibranch tree
	//	m_multiBody->addChild(m_multiBody2);
	//	m_multiBody2->setParentCustomBody(m_multiBody);

	//	btQuaternion baseQ = btQuaternion(0, 0, 0, 1);
	//	pMultiBody->setWorldToBaseRot(baseQ);
	//	if (m_multiBody2->getParentCustomBody() != nullptr)
	//	{
	//		customMultiBody* parent = m_multiBody2->getParentCustomBody();
	//		pMultiBody->setBasePos(parent->m_tailPos);
	//	}
	//	else
	//	{
	//		pMultiBody->setBasePos(btVector3(0, 0, 0));
	//	}

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
	//		if (i == 0 && m_multiBody2->getParentCustomBody() != nullptr)
	//		{
	//			btQuaternion q = btQuaternion(btVector3(0, 1, 0).normalized(), 30 * SIMD_PI / 180.f);
	//			//btQuaternion q = btQuaternion(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
	//			//btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			btQuaternion tail = m_multiBody2->getParentCustomBody()->m_tailRot;
	//			q = tail * q;
	//			pMultiBody->setJointPosMultiDof(i, q);
	//		}
	//		else
	//		{
	//			btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			pMultiBody->setJointPosMultiDof(i, q);
	//		}
	//	}
	//	//for (int i = 0; i < pMultiBody->getNumLinks(); i++)
	//	//{
	//	//	btQuaternion q(btVector3(1, 0, 0).normalized(), 20 * SIMD_PI / 180.f);
	//	//	pMultiBody->setJointPosMultiDof(i, q);
	//	//}

	//	// init balance
	//	pMultiBody->m_balanceRot.resize(pMultiBody->getNumLinks());
	//	for (int i = 0; i < numLink; i++)
	//	{
	//		if (i == 0 && m_multiBody2->getParentCustomBody() != nullptr)
	//		{
	//			btQuaternion q = btQuaternion(btVector3(0, 1, 0).normalized(), 30 * SIMD_PI / 180.f);
	//			pMultiBody->m_balanceRot[i] = q;
	//		}
	//		else
	//		{
	//			btQuaternion q = btQuaternion(0, 0, 0, 1);
	//			pMultiBody->m_balanceRot[i] = q;
	//		}
	//	}

	//	// use m_ks from main skeleton

	//	addColliders(pMultiBody, m_dynamicsWorld, joint_lengths);
	//}

	//// init p2p constraint
	//{
	//	btVector3 pointInA = m_multiBody->worldPosToLocal(m_numLinks - 1, m_multiBody->m_tailPos);
	//	btVector3 pointInB = m_multiBody2->worldPosToLocal(0, m_multiBody->m_tailPos);
	//	btMultiBodyPoint2Point* p2p = new btMultiBodyPoint2Point(m_multiBody, m_numLinks - 1, m_multiBody2, 0, pointInA, pointInB);
	//	p2p->setMaxAppliedImpulse(100);
	//	p2p->setErp(0.8);
	//	m_dynamicsWorld->addMultiBodyConstraint(p2p);
	//}


	/////////////////////////////////////////////////////////////////
	// construct the tree
	/////////////////////////////////////////////////////////////////
    m_multiBody->constructTreeList();

    /////////////////////////////////////////////////////////////////
    // construct the marker
    /////////////////////////////////////////////////////////////////
  //  {
		//btCollisionShape* colShape = new btSphereShape(btScalar(0.02 * LENGTH_RATIO));

  //      /// Create Dynamic Objects
  //      btTransform startTransform;
  //      startTransform.setIdentity();
		//startTransform.setOrigin(m_multiBody->m_tailPos);

  //      btScalar mass(0.f);

  //      //rigidbody is dynamic if and only if mass is non zero, otherwise static
  //      bool isDynamic = (mass != 0.f);

  //      btVector3 localInertia(0, 0, 0);
  //      if (isDynamic)
  //          colShape->calculateLocalInertia(mass, localInertia);

  //      //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
  //      btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
  //      btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
  //      m_marker = new btRigidBody(rbInfo);
  //      m_marker->setRestitution(0.0);
  //      m_dynamicsWorld->addRigidBody(m_marker, NOTHING, NOTHING);

  //      if (!WIRE_FRAME)
  //      {
  //          m_guiHelper->createCollisionShapeGraphicsObject(colShape);
  //          btVector4 color(1, 0, 0, 1);
  //          m_guiHelper->createCollisionObjectGraphicsObject(dynamic_cast<btCollisionObject*>(m_marker), color);
  //      }
  //  }

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

void Skeleton::applySpringForce(customMultiBody* mb)
{
    // if reduce_factor > 0.8, it looks unnatural. 
	btScalar reduce_factor = 0.6;

    btAlignedObjectArray<btVector3> fbXaxisTorques;
	fbXaxisTorques.resize(mb->getNumLinks(), btVector3(0, 0, 0));

	btAlignedObjectArray<btVector3> XaxisTorqueDirs;
	XaxisTorqueDirs.resize(mb->getNumLinks(), btVector3(0, 0, 0));

    btAlignedObjectArray<btVector3> fbRadialTorques;
	fbRadialTorques.resize(mb->getNumLinks(), btVector3(0, 0, 0));

    btAlignedObjectArray<btVector3> RadialTorqueDirs;
	RadialTorqueDirs.resize(mb->getNumLinks(), btVector3(0, 0, 0));

    const btScalar max_torque = 0.1 * MASS_RATIO;
    const btVector3 joint_dir = btVector3(1, 0, 0);

    btQuaternion tailRot(0, 0, 0, 1);
	btQuaternion tailBalanceQ(0, 0, 0, 1);
	if (mb->getParentCustomBody() != nullptr)
	{
		tailRot = mb->getParentCustomBody()->m_tailRot;		
		tailBalanceQ = mb->getParentCustomBody()->m_tailBalanceRot;
    }
	btQuaternion prevXQ = tailRot;
	btQuaternion prevQ = tailBalanceQ;

    // forwards
	btVector3 last_c(0, 0, 0);
	
    for (int i = 0; i < mb->getNumLinks(); ++i)
    {
		btVector3 c = mb->getLink(i).m_cachedWorldTransform.getOrigin();

        if (i == (mb->getNumLinks() - 1))
			last_c = c;

        // cause joint_dir is a unit vector, so after localDirToWorld it does not need normalization.
		btVector3 dir = mb->localDirToWorld(i, joint_dir);
        btVector3 origin = c + (-dir) * joint_lengths[i] * 0.5;
		btQuaternion balanceQ = mb->m_balanceRot[i];
		btVector3 torque(0, 0, 0);

        // x axis
        {
            btQuaternion currentQuat(mb->getJointPosMultiDof(i)[0], mb->getJointPosMultiDof(i)[1], mb->getJointPosMultiDof(i)[2],
							mb->getJointPosMultiDof(i)[3]);
            btVector3 angleDiff;
			// TODO remove i == 0
            btQuaternion cq = (i==0 ? btQuaternion(0,0,0,1) : prevXQ) * currentQuat;
			btQuaternion tq = prevXQ * balanceQ;
			btQuaternion middle = cq.slerp(tq, 0.5);
			btQuaternion relRot = middle.inverse() * tq;
			btGeneric6DofSpring2Constraint::matrixToEulerXYZ(btMatrix3x3(relRot), angleDiff);
			btScalar angle = angleDiff.x();

            if (fabs(angle) > 1 * SIMD_PI / 180.f)
			{
				XaxisTorqueDirs[i] = -dir;
				//dampingFunc* pDampFunc = static_cast<dampingFunc*>(mb->getLink(i).m_collider->getUserPointer());
				//btScalar strength = pDampFunc->getChange(angle);
				torque = XaxisTorqueDirs[i] * m_Ks[i] * (joint_lengths[i] * 0.5) * angle;
				torque *= 5;

				if (torque.safeNorm() > max_torque)
				{
					torque *= max_torque / torque.safeNorm();
				}
				mb->addLinkTorque(i, torque);
				fbXaxisTorques[i] += -torque;
			}

            if (WIRE_FRAME)
			{
				/*if (!torque.fuzzyZero())
				    m_guiHelper->getRenderInterface()->drawLine(c, c + torque.safeNormalize() * 0.05, btVector4(1, 0, 0, 1), btScalar(2));*/
			}
			prevXQ = (i == 0 ? btQuaternion(0, 0, 0, 1) : prevXQ) * currentQuat;
		}

        // radial axis
        {
			btVector3 balanceP = quatRotate(prevQ * balanceQ, joint_dir * joint_lengths[i]) * 0.5 + origin;
			btScalar angle = btAcos((balanceP - origin).normalized().dot(dir));
			if (angle > SIMD_PI / 6)
			{
				angle = SIMD_PI / 6;
			}
			btScalar weight = 0.85;
			angle *= weight;
			//if (angle > 1 * SIMD_PI / 180.f)
			//{
			
			btVector3 rtdir = dir.cross((balanceP - origin));
			if (!rtdir.fuzzyZero())
			{
				rtdir.normalize();
				//if (rtdir.z() < 0)
				//{
				//	printf("\n");
				//}
				if (!mb->m_firstspringcalc)
				{
					if (btIsNegative(mb->m_rtdirs[i].dot(rtdir)))
					{
						mb->m_signs[i] = -mb->m_signs[i];
					}
					if (rtdir.length() > 1e-7)
					{
						mb->m_rtdirs[i] = rtdir;
					}
					angle *= mb->m_signs[i] * angle;
				}
				RadialTorqueDirs[i] = rtdir;

				dampingFunc* pDampFunc = static_cast<dampingFunc*>(mb->getLink(i).m_collider->getUserPointer());
				btScalar strength = pDampFunc->getChange(angle);
				int brake = 1;
				if (btIsNegative(strength * angle))
				{
					brake = -1;
				}
				//btScalar strength = 1.96;
				printf("joint %d: angle %f, strength %f\n", i, angle, strength);

				btScalar ml = 0.0;
				for (int j = i; j < mb->getNumLinks(); j++)
				{
					ml += (j - i) + 0.5;
				}
				torque = RadialTorqueDirs[i] * m_Ks[i] * brake * fabs(strength) * sin(fabs(angle)) * (joint_lengths[i] * ml) * 20;
				/*if (brake < 0)
				{
					printf("joint %d: brake\n", i);
					torque = btVector3(0, 0, 0);
				}
				else
				{
					printf("joint %d: continue\n", i);
				}*/

				// limit the max torque
				if (torque.safeNorm() > max_torque)
				{
					torque *= max_torque / torque.safeNorm();
				}
				mb->addLinkTorque(i, torque);

				//btVector3 rf = rtdir.cross(dir) * power * sin(angle);
				//mb->addLinkForce(i, rf);
				//btVector3 lf = dir * power * cos(angle);
				//mb->addLinkForce(i, lf);

				fbRadialTorques[i] += -torque;

				if (WIRE_FRAME)
				{
					if (!torque.fuzzyZero())
					{
						//m_guiHelper->getRenderInterface()->drawLine(c, c + lf.safeNormalize() * 0.05, btVector4(0, 1, 0, 1), btScalar(2));
						m_guiHelper->getRenderInterface()->drawPoint(balanceP, btVector4(0, 0, 1, 1), btScalar(5));
					}
				}
			}
			//}
			
			prevQ = prevQ * balanceQ;
		}
    }

	if (mb->m_firstspringcalc)
	{
		mb->m_rtdirs.resize(mb->getNumLinks());
		mb->m_signs.resize(mb->getNumLinks());
		for (size_t i = 0; i < mb->getNumLinks(); i++)
		{
			mb->m_rtdirs[i] = RadialTorqueDirs[i];
			mb->m_signs[i] = 1;
		}
		mb->m_firstspringcalc = false;
	}

    // backwards
	std::vector<customMultiBody*> childs = mb->getChilds();
	// x axis feedback
	btVector3 tailXaxisTorque(0, 0, 0);
	for (int i = 0; i < childs.size(); i++)
	{
		tailXaxisTorque += childs[i]->m_fbXaxisTorque;
	}
	for (int i = mb->getNumLinks() - 1; i >= 0; --i)
	{
		if (i == mb->getNumLinks() - 1)
		{
			btVector3 fb = tailXaxisTorque * reduce_factor;
			fb = fb.dot(XaxisTorqueDirs[i]) * XaxisTorqueDirs[i];
			//mb->addLinkTorque(i, fb);
		}
		else
		{
			btVector3 fb = fbXaxisTorques[i + 1] * reduce_factor;
			fb = fb.dot(XaxisTorqueDirs[i]) * XaxisTorqueDirs[i];
			//mb->addLinkTorque(i, fb);
		}
	}

    // radial feedback
	btVector3 tailRadialTorque(0, 0, 0);
	for (int i = 0; i < childs.size(); i++)
	{
		tailRadialTorque += childs[i]->m_fbRadialTorque;
	}

    for (int i = mb->getNumLinks() - 1; i >= 0; --i)
	{
		if (i == mb->getNumLinks() - 1)
		{
			btVector3 fb = tailRadialTorque * reduce_factor;
			fb = fb.dot(RadialTorqueDirs[i]) * RadialTorqueDirs[i];
			//mb->addLinkTorque(i, fb);
		}
		else
		{
			btVector3 fb = fbRadialTorques[i + 1] * reduce_factor;
			fb = fb.dot(RadialTorqueDirs[i]) * RadialTorqueDirs[i];
			//mb->addLinkTorque(i, fb);
        }
	}
    
    mb->m_fbRadialTorque = fbRadialTorques[0];
	mb->m_fbXaxisTorque = fbXaxisTorques[0];

    if (WIRE_FRAME)
	{
		if (!tailRadialTorque.fuzzyZero())
		{
			//m_guiHelper->getRenderInterface()->drawLine(last_c, last_c + tailRadialTorque.safeNormalize() * 0.05, btVector4(0, 0, 1, 1), btScalar(2));
		}
	}
}

void Skeleton::addColliders(customMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const float* joint_lengths)
{
    btAlignedObjectArray<btQuaternion> world_to_local;
    world_to_local.resize(pMultiBody->getNumLinks() + 1);

    btAlignedObjectArray<btVector3> local_origin;
    local_origin.resize(pMultiBody->getNumLinks() + 1);
    world_to_local[0] = btQuaternion(0,0,0,1);
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

        if (i == (pMultiBody->getNumLinks() -1))
		{
			pMultiBody->m_tailRot = btQuaternion(quat[0], quat[1], quat[2], quat[3]);
			pMultiBody->m_tailPos = posr + quatRotate(pMultiBody->m_tailRot, btVector3(joint_lengths[i] / 2.0f, 0, 0));
        }

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

		//pdDamp* damp = new pdDamp(5.0);
		//damp->setPID(0.5, 0.2, 0.7);
		//damp->setMaxIOutput(SIMD_PI * 2.0);
		//damp->setOutputLimits(-SIMD_PI * 1.0, SIMD_PI * 1.0);
		//damp->setMaxError(SIMD_PI * 1.0);
		//damp->setOutputFilter(0.2);

		//linearDamp* damp = new linearDamp(5.0);

		pdDamp* damp = new pdDamp(1.0, 3.14, 0.4, 30, 0.2, 0.2);

		col->setUserPointer(static_cast<void*>(damp));

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

    demo->m_multiBody->updateTailRots();
	for (int i = demo->m_multiBody->m_tree.size() - 1; i >= 0; i--)
	{
		customMultiBody* curr = demo->m_multiBody->m_tree[i];
		demo->applySpringForce(curr);
    }

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

    //printf("step: %d\n", m_step);

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
        //const char* gPngFileName = "multibody";
        //sprintf(mfileName, "%s_%d.png", gPngFileName, m_step);
        //this->m_guiHelper->getAppInterface()->dumpNextFrameToPng(mfileName);
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
	printf("step %d\n", m_step);
    if (m_step == 28)
    {
        printf("step %d\n", m_step);
    }
}


class CommonExampleInterface* Dof6Spring2CreateFunc(CommonExampleOptions& options)
{
//    return new Dof6Spring2Setup(options.m_guiHelper);
	return new Skeleton(options.m_guiHelper);
}
