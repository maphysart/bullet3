#include <math.h>
#include <limits>
#include <deque>

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

#define MAX_SIZE 3

#define ACC_ARRAY_SIZE 20

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

    btVector3 getPos(float time, float deltaTime){
        btVector3 pos = btVector3(0, 0, cosOffset(amp, T, phase, time) * exp(-time / (4* T)));
        return pos;
    }

    btScalar amp;
    btScalar T;
    btScalar phase;
};


struct CircleOscilator
{
    CircleOscilator()
    {
        T = 30.0f;
        omega = SIMD_2_PI / T;
        angle = 0.0f;
        clockwise = true;
        radius = 2;
    }

    btVector3 getPos(float time, float deltaTime){
        if ( angle > (SIMD_PI / 4.0) )
        {
            clockwise = false;
        }
        if ( angle < (-SIMD_PI / 4.0) )
        {
            clockwise = true;
        }
        angle += (clockwise ? 1 : -1) * omega * deltaTime * exp(-time / (2 * T));
//        angle += (clockwise ? 1 : -1) * omega * deltaTime;
        return btVector3(radius * sin(angle), 0, radius * cos(angle));
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
    btRigidBody* m_collider;

    int m_solverType;
    btScalar m_time;
    int m_step;
    int m_numLinks;
    btAlignedObjectArray<btQuaternion> m_balanceRot;
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

public:
    Skeleton(struct GUIHelperInterface* helper);
    virtual ~Skeleton();
    virtual void initPhysics();
    virtual void stepSimulation(float deltaTime);
    virtual void resetCamera()
    {
        float dist = 10;
        float pitch = -30;
        float yaw = 90;
        float targetPos[3] = {0, 0, 0};
        m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
    }
    void addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents);
    void moveCollider(const btVector3& pos);
    void applySpringForce(float time, float deltaTime);
    void applyBaseLinearDragForce(const btVector3& dir);
    void applyBaseCentrifugalForce(float deltaTime, int m_step, const btVector3& pos);
    static void OnInternalTickCallback(btDynamicsWorld* world, btScalar timeStep);
    void applyGravityForce(float deltaTime);
    void getLinearAcc(float deltaTime, int m_step, const btVector3& pos);
    void limitMaxTwist(float max_angle);
};

Skeleton::Skeleton(struct GUIHelperInterface* helper)
        : CommonMultiBodyBase(helper)
{
    m_time = btScalar(0.0);
    m_step = 0;
    m_numLinks = 16;
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

    btVector3 init_pos = m_move.getPos(0, 0);

    /////////////////////////////////////////////////////////////////
    // construct the skeleton
    /////////////////////////////////////////////////////////////////
    btVector3 linkHalfExtents(0.1, 0.5, 0.1);
    btVector3 baseHalfExtents(0.0, 0.0, 0.0);

    btVector3 baseInertiaDiag(0.f, 0.f, 0.f);
    float baseMass = 1.f;

    if (baseMass)
    {
        btCollisionShape* pTempBox = new btSphereShape(btScalar(0.1));
        pTempBox->calculateLocalInertia(baseMass, baseInertiaDiag);
        delete pTempBox;
    }

    customMultiBody* pMultiBody = new customMultiBody(m_numLinks, baseMass, baseInertiaDiag, !floating, canSleep);
    //pMultiBody->useRK4Integration(true);

    // set base position
    btQuaternion baseQ(0.f, 0.f, 0.f, 1.f);
    pMultiBody->setWorldToBaseRot(baseQ);
    btVector3 basePos = btVector3(0.0, 0.0, m_radius);
    pMultiBody->setBasePos(init_pos);

    m_multiBody = pMultiBody;

    //y-axis assumed up
    btVector3 currentPivotToCurrentCom(0, -linkHalfExtents[1], 0);
    for (int i = 0; i < m_numLinks; ++i)
    {
        float linkMass = 1.f;
        btVector3 linkInertiaDiag(0.f, 0.f, 0.f);
        btCollisionShape* shape = 0;
        {
            shape = new btBoxShape(linkHalfExtents);
        }
        shape->calculateLocalInertia(linkMass, linkInertiaDiag);
        delete shape;

        btVector3 parentComToCurrentCom;
        if (i == 0){
            parentComToCurrentCom = btVector3(0, -linkHalfExtents[1] * 1, 0);
        }
        else{
            parentComToCurrentCom = btVector3(0, -linkHalfExtents[1] * 2, 0);
        }
        pMultiBody->setupSpherical(i, linkMass, linkInertiaDiag, i - 1,
                btQuaternion(0.f, 0.f, 0.f, 1.f),
                                  (parentComToCurrentCom - currentPivotToCurrentCom),
                                  currentPivotToCurrentCom, !selfCollide);
    }

    // init params
    pMultiBody->finalizeMultiDof();
    m_dynamicsWorld->addMultiBody(pMultiBody);
    pMultiBody->setCanSleep(canSleep);
    pMultiBody->setHasSelfCollision(selfCollide);
    pMultiBody->setUseGyroTerm(gyro);
    if (damping)
    {
        // TODO set linear and angular damp for each joint
        pMultiBody->setLinearDamping(0.1f);
        pMultiBody->setAngularDamping(0.7f);
    }

    // init pose
    btScalar angle_damp = 0.0;
    btScalar angle;

    angle = 0 * SIMD_PI / 180.f;
    for ( int i = 0; i < pMultiBody->getNumLinks(); i++ )
    {
        btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
        pMultiBody->setJointPosMultiDof(i, q);
        angle *= angle_damp;
    }

    // init balance
    m_balanceRot.resize(pMultiBody->getNumLinks());
    angle = 0 * SIMD_PI / 180.f;
    for (int i = 0; i < m_numLinks; i++)
    {
        if ( i == 0)
        {
            angle = 45 * SIMD_PI / 180.f;
            btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
            m_balanceRot[i] = q;
        } else{
            btQuaternion q(btVector3(1, 0, 0).normalized(), 0);
            m_balanceRot[i] = q;
        }
    }

    m_Ks.resize(pMultiBody->getNumLinks());
    float ks_damp = 0.2;
    float ks_init = 0.2;
    for (int i = 0; i < m_numLinks; i++)
    {
        m_Ks[i] = ks_init / pow(i+1, 1.0);
    }

    addColliders(pMultiBody, m_dynamicsWorld, baseHalfExtents, linkHalfExtents);


    /////////////////////////////////////////////////////////////////
    // construct the box
    /////////////////////////////////////////////////////////////////
    {
        btCollisionShape* colShape = new btSphereShape(btScalar(0.1));

        /// Create Dynamic Objects
        btTransform startTransform;
        startTransform.setIdentity();
        startTransform.setOrigin(btVector3(0.0, 0.0, 0.0));

        btScalar mass(0.f);

        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);

        btVector3 localInertia(0, 0, 0);
        if (isDynamic)
            colShape->calculateLocalInertia(mass, localInertia);

        //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
        btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
        m_collider = new btRigidBody(rbInfo);
        m_collider->setRestitution(0.0);
        m_dynamicsWorld->addRigidBody(m_collider, BONE_BODY, NOTHING);
        int id = m_collider->getCompanionId();

        if (!WIRE_FRAME)
        {
            m_guiHelper->createCollisionShapeGraphicsObject(colShape);
            btVector4 color(1, 0, 0, 1);
            m_guiHelper->createCollisionObjectGraphicsObject(dynamic_cast<btCollisionObject*>(m_collider), color);
        }
    }

    // init p2p constraint
    {
        btVector3 pointInA = pMultiBody->worldPosToLocal(0, init_pos);
        btVector3 pointInB = init_pos;
        btMatrix3x3 frameInA;
        btMatrix3x3 frameInB;
        frameInA.setIdentity();
        frameInB.setIdentity();
        m_p2p = new btMultiBodyPoint2Point(pMultiBody, 0, 0, pointInA, pointInB);
        m_p2p->setMaxAppliedImpulse(10);
        m_p2p->setErp(0.8);
        m_dynamicsWorld->addMultiBodyConstraint(m_p2p);
    }

    // init motor constraint
//    {
//        btScalar angle = 45 * SIMD_PI / 180.f;
//        float damp = 0.0;
//        for (int i = 0; i < pMultiBody->getNumLinks(); i++) {
//            btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
//            btMultiBodySphericalJointMotor *motor = new btMultiBodySphericalJointMotor(pMultiBody, i, btScalar(10));
//            motor->setPositionTarget(q, 0.1);
//            motor->setErp(btScalar(0.1));
//            motor->finalizeMultiDof();
//            angle *= damp;
//            motors.push_back(motor);
//        }
//    }

    m_dynamicsWorld->setGravity(btVector3(0, -0.0, 0));

    btContactSolverInfo& si = m_dynamicsWorld->getSolverInfo();
    si.m_numIterations = 10;
    si.m_globalCfm = 0.05f;
    si.m_erp2 = 0.1f;

//    m_multiBody->setMaxCoordinateVelocity(2.0);
    m_multiBody->setMaxOmega(1.0f);
}

void Skeleton::addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents)
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
    btTransform tr;
    tr.setIdentity();
    tr.setOrigin(pos);
    tr.setRotation(btQuaternion(0.0, 0.0, 0.0, 1.0));
    m_collider->setWorldTransform(tr);
}

void Skeleton::applySpringForce(float time, float deltaTime){
    float max_force = 0.2;
    std::vector<btScalar> temp;
//    float max_time = 0.5;
//    temp = m_Ks;

    for (int i = 0; i < m_multiBody->getNumLinks(); ++i) {
        btQuaternion currentQuat(m_multiBody->getJointPosMultiDof(i)[0],
                                 m_multiBody->getJointPosMultiDof(i)[1],
                                 m_multiBody->getJointPosMultiDof(i)[2],
                                 m_multiBody->getJointPosMultiDof(i)[3]);
        btQuaternion relRot = currentQuat * m_balanceRot[i].inverse();
        btVector3 angleDiff;
        btGeneric6DofSpring2Constraint::matrixToEulerXYZ(btMatrix3x3(relRot), angleDiff);
        btVector3 spring;
        for (int d = 0; d < m_multiBody->getLink(i).m_dofCount; d++) {
            spring[d] = m_Ks[i] * angleDiff[d];
        }
        if ( spring.norm() > max_force )
        {
            spring = spring.normalized() * max_force;
        }
//        m_multiBody->addLinkForce(i, btVector3(0, 0, -1) * m_Ks[i]);
        m_multiBody->addLinkTorque(i, spring);
//        spring[1] = 0.0f;
//        m_multiBody->addJointTorqueMultiDof(i, spring);
//        m_multiBody->addJointTorqueMultiDof(i, &spring[0]);
//        for (int d = 0; d < m_multiBody->getLink(i).m_dofCount; d++) {
//            m_multiBody->addJointTorqueMultiDof(i, d, spring[d]);
//        }
    }
}

void Skeleton::applyGravityForce(float deltaTime)
{
    float damp = 1.0; // if damp is set to 0.2 or more, the tail will wind up in horizontal oscillating motion
    btVector3 g = btVector3(0, -0.2, 0);
    for (int i = 0; i < m_numLinks; ++i) {
        m_multiBody->addLinkForce(i, g * m_multiBody->getLink(i).m_mass / pow(i+1, damp));
    }
}

void Skeleton::getLinearAcc(float deltaTime, int m_step, const btVector3& pos)
{
    if ( positions.size() == MAX_SIZE)
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
    printf("acc %f\n", m_currAcc.norm());
    printf("vel %f\n", m_currVel.norm());
}

void Skeleton::applyBaseLinearDragForce(const btVector3& dir)
{
    btVector3 currBaseVel;
    float scaling = 5;
    float damp = 1.0; // if damp is set to 0.2 or more, the tail will wind up in horizontal oscillating motion

    // TODO replace to avg acc
    if ( m_avgAcc > 1e-2 ) {
        btVector3 g = dir * m_avgAcc * scaling;
        for (int i = 0; i < m_numLinks; ++i) {
            m_multiBody->addLinkForce(i, g * m_multiBody->getLink(i).m_mass / pow(i+1, damp));
        }
    }
}

void Skeleton::applyBaseCentrifugalForce(float deltaTime, int m_step, const btVector3& pos)
{
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

//    // calculate and apply the impulse, the damp use the difference between prev pos and curr pos
//    applySpringForce(m_time, deltaTime);
    printf("step: %d\n", m_step);

//    btVector3 basePos = m_move.getPos(m_time, deltaTime);
    btVector3 basePos = m_move.getPos(m_time, deltaTime);

    getLinearAcc(deltaTime, m_step, basePos);

    if ( accs.size() == ACC_ARRAY_SIZE ) {
        float sum = 0.0f;
        for ( int i = 0; i < accs.size(); i++ )
            sum += accs[i].norm();
        m_avgAcc = sum / ACC_ARRAY_SIZE;
    }

    /*
     * try to combine two constraint methods, when base is not accelerated, we use spring force,
     * if the base acc is big, use btMultiBodySphericalJointMotor, but the problem is that
     * the kp setPositionTarget(q, kp) in btMultiBodySphericalJointMotor does not work even set
     * to 0.
     */
//    if ( accs.size() == ACC_ARRAY_SIZE )
//    {
//        float avg_acc = 0.0f;
//        // if average acc in several frames is less than thresh hold, use drag force to simulate
//        for ( int i = 0; i < accs.size(); i++ )
//            avg_acc += accs[i].norm();
//        avg_acc /= ACC_ARRAY_SIZE;
//        if ( avg_acc <= 0.01 )
//        {
//            if (use_constraint)
//            {
//                for ( int i = 0; i < motors.size(); i++ )
//                {
//                    m_dynamicsWorld->removeMultiBodyConstraint(motors[i]);
//                    delete motors[i];
//                }
//                motors.clear();
//                use_constraint = false;
//            }
//            applySpringForce(m_time, deltaTime);
//        }
//        else {
//            {
//                if ( !use_constraint ){
//                    btScalar angle = 45 * SIMD_PI / 180.f;
//                    float damp = 0.0;
//                    for (int i = 0; i < m_multiBody->getNumLinks(); i++) {
//                        btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
//                        btMultiBodySphericalJointMotor *motor = new btMultiBodySphericalJointMotor(m_multiBody, i, btScalar(0.1));
//                        motor->setPositionTarget(q, 0.0);
//                        motor->setErp(btScalar(0.1));
//                        m_multiBody->getLink(i).m_userPtr = motor;
//                        m_dynamicsWorld->addMultiBodyConstraint(motor);
//                        motor->finalizeMultiDof();
//                        angle *= damp;
//                        motors.push_back(motor);
//                    }
//                    use_constraint = true;
//                }
//            }
//        }
//    }


    btVector3 dir = btVector3(0, -1, 0);

    applyBaseLinearDragForce(dir);

//    applyBaseCentrifugalForce(deltaTime, m_step, basePos);

    applySpringForce(m_time, deltaTime);

    // p2p
    m_p2p->setPivotInB(basePos);

    moveCollider(basePos);

    // capture the frames
//    {
//        const char* gPngFileName = "multibody";
//        sprintf(mfileName, "%s_%d.png", gPngFileName, m_step);
//        this->m_guiHelper->getAppInterface()->dumpNextFrameToPng(mfileName);
//    }

    // step and update positions
    int substeps = 2;
    btScalar  fixedTimeStep = deltaTime / substeps;
    m_dynamicsWorld->stepSimulation(deltaTime, substeps, fixedTimeStep);
    if ( WIRE_FRAME ) {
        m_dynamicsWorld->debugDrawWorld();
    }

    btScalar max_angle = 70 * SIMD_PI / 180.f;
    limitMaxTwist(max_angle);

    m_time += deltaTime;
    m_step += 1;

    if ( m_step == 121 )
    {
        printf("step\n");
    }
}


class CommonExampleInterface* Dof6Spring2CreateFunc(CommonExampleOptions& options)
{
//    return new Dof6Spring2Setup(options.m_guiHelper);
	return new Skeleton(options.m_guiHelper);
}
