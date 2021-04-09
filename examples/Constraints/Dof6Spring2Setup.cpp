#include <math.h>
#include <limits>

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

struct lineOscilator
{
    lineOscilator()
    {
        amp = 0.5f;
        T = 50.0f;
        phase =  0;
    }

    btVector3 getPos(float time, float deltaTime){
        btVector3 pos = btVector3(0, 0, cosOffset(amp, T, phase, time));
        return pos;
    }

    btScalar amp;
    btScalar T;
    btScalar phase;
};

struct Skeleton : public CommonMultiBodyBase
{
    btMultiBody* m_multiBody;
    btRigidBody* m_collider;

    int m_solverType;
    btScalar m_time;
    int m_step;
    int m_numLinks;
    btVector3 m_prevBaseVel;
    btVector3 m_prevBasePos;
    btAlignedObjectArray<btQuaternion> m_balanceRot;
    btAlignedObjectArray<btScalar> m_Ks;
    gravityGenerator m_g;
    btScalar m_linearDragEffect;
    btScalar m_centrifugalDragEffect;
    float m_angle;
    bool m_clockwise;
    btMultiBodyPoint2Point* m_p2p;
    btScalar m_radius;
    lineOscilator m_move;

public:
    Skeleton(struct GUIHelperInterface* helper);
    virtual ~Skeleton();
    virtual void initPhysics();
    virtual void stepSimulation(float deltaTime);
    virtual void resetCamera()
    {
        float dist = 8;
        float pitch = -21;
        float yaw = 90;
        float targetPos[3] = {0, 0, 0};
        m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
    }
    void addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents);
    btTransform transformBase(btScalar time, btScalar deltaTime);
    void moveCollider(const btVector3& pos);
    void applySpringForce(float deltaTime);
    void applyBaseLinearDragForce(float deltaTime, int m_step, const btVector3& pos);
    void applyBaseCentrifugalForce(float deltaTime, int m_step, const btVector3& pos);
    static void OnInternalTickCallback(btDynamicsWorld* world, btScalar timeStep);
    void applyGravityForce(float deltaTime);
};

Skeleton::Skeleton(struct GUIHelperInterface* helper)
        : CommonMultiBodyBase(helper)
{
    m_time = btScalar(0.0);
    m_step = 0;
    m_numLinks = 8;
    m_solverType = 0;
    m_prevBaseVel = btVector3(0.0, 0.0, 0.0);
    m_prevBasePos = btVector3(0, 0,0);
    m_g.m_gravity = -0.08;
    m_linearDragEffect = btScalar(25.0);
    m_centrifugalDragEffect = btScalar(0.3);
    m_angle = 0.0;
    m_clockwise = true;
    m_radius = 2.0;
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
    const bool damping = false;   // disable bullet internal damp
    const bool gyro = false;
    const bool canSleep = false;
    const bool selfCollide = false;

    m_move.amp = 2;
    m_move.T = 10;
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

    btMultiBody* pMultiBody = new btMultiBody(m_numLinks, baseMass, baseInertiaDiag, !floating, canSleep);
    //pMultiBody->useRK4Integration(true);

    // set base position
    btQuaternion baseQ(0.f, 0.f, 0.f, 1.f);
    pMultiBody->setWorldToBaseRot(baseQ);
    btVector3 basePos = btVector3(0.0, 0.0, m_radius);
    pMultiBody->setBasePos(init_pos);

    m_prevBasePos = init_pos;
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
        pMultiBody->setLinearDamping(0.05f);
        pMultiBody->setAngularDamping(0.1f);
    }

    // init pose
    btScalar angle_damp = 0.0;
    btScalar angle = 0 * SIMD_PI / 180.f;
    std::vector<btQuaternion> qs;
    for ( int i = 0; i < pMultiBody->getNumLinks(); i++ )
    {
        btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
        qs.push_back(q);
        pMultiBody->setJointPosMultiDof(i, q);
        angle *= angle_damp;
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
        m_dynamicsWorld->addMultiBodyConstraint(m_p2p);
    }

    // init motor constraint
    {
        btScalar angle = 0 * SIMD_PI / 180.f;
        float damp = 0.0;
        for (int i = 0; i < pMultiBody->getNumLinks(); i++) {
            btQuaternion q(btVector3(1, 0, 0).normalized(), angle);
            btMultiBodySphericalJointMotor *motor = new btMultiBodySphericalJointMotor(pMultiBody, i, btScalar(10));
            motor->setPositionTarget(q, 0.1);
            motor->setErp(btScalar(0.1));
            pMultiBody->getLink(i).m_userPtr = motor;
            m_dynamicsWorld->addMultiBodyConstraint(motor);
            motor->finalizeMultiDof();

            angle *= damp;
        }
    }

    m_dynamicsWorld->setGravity(btVector3(0, -0.0, 0));
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

btTransform Skeleton::transformBase(btScalar time, btScalar deltaTime){
    btTransform trans;
    trans.setIdentity();

    // btMultiBodyFixedConstraint does not work as we like, cause it is not totally fixed to the body.
    // move the base of the chain does not apply torque on the multibody chains
//    btScalar amp = 0.5f;
//    btScalar T = 50.0f;
//    btScalar phase =  SIMD_PI / 2.0;
//    btVector3 basePos = btVector3(cosOffset(amp, T, phase, time), 0, 0);
//    trans.setOrigin(basePos);

    btScalar cycle = 100.0f;
    btScalar omega = SIMD_2_PI / cycle;
    if ( m_angle > (SIMD_PI / 4.0) )
    {
        m_clockwise = false;
    }
    if ( m_angle < (-SIMD_PI / 4.0) )
    {
        m_clockwise = true;
    }
    omega *= m_clockwise ? 1 : -1;
    m_angle += omega * deltaTime;
    btQuaternion q;
    q.setRotation(btVector3(0, 1, 0), m_angle);
    trans.setRotation(q);

    btScalar radius = 1.0;
    btVector3 basePos = btVector3(m_radius * sin(m_angle), 0, m_radius * cos(m_angle));
    trans.setOrigin(basePos);

    m_multiBody->setBaseWorldTransform(trans);

    return trans;
}

void Skeleton::moveCollider(const btVector3& pos){
    btTransform tr;
    tr.setIdentity();
    tr.setOrigin(pos);
    tr.setRotation(btQuaternion(0.0, 0.0, 0.0, 1.0));
    m_collider->setWorldTransform(tr);
}

void Skeleton::applySpringForce(float deltaTime){
    for (int i = 0; i < m_multiBody->getNumLinks(); ++i) {
        if (m_Ks[i] > 0.0)
        {
            btScalar* q = m_multiBody->getJointPosMultiDof(i);
            btQuaternion curr_q = btQuaternion(q[0], q[1], q[2], q[3]);
            btScalar _ks = m_Ks[i];
            // TODO ?
//            _ks = adjustKs(m_Ks[i], m_multiBody->getLink(i).m_mass, deltaTime);
            btVector3 spring(m_balanceRot[i][0] - curr_q[0],m_balanceRot[i][1] - curr_q[1], m_balanceRot[i][2] - curr_q[2]);
            spring *= _ks;
            if (spring.length() > MAX_SPRING_FORCE) {
                printf("spring force %f is too large\n", spring.length());
                spring = spring.normalized() * MAX_SPRING_FORCE;
            }
            // apply the torque
            for (int d = 0; d < m_multiBody->getLink(i).m_dofCount; d++) {
                m_multiBody->addJointTorqueMultiDof(i, d, spring[d]);
            }
        }
    }
}

void Skeleton::applyGravityForce(float deltaTime)
{
    btVector3 g = btVector3(0, -10, 0);
    for (int i = 0; i < m_numLinks; ++i) {
        btVector3 orient = m_multiBody->getLink(i).m_dVector;
        orient = m_multiBody->localDirToWorld(i, orient);
        btVector3 torque = btCross(orient, g * m_multiBody->getLink(i).m_mass);
//        torque *= deltaTime;
        m_multiBody->addLinkTorque(i, torque);
    }
//    for (int i = 0; i < m_numLinks; ++i) {
//        m_multiBody->addLinkForce(i, g * m_multiBody->getLink(i).m_mass);
//    }
}

void Skeleton::applyBaseLinearDragForce(float deltaTime, int m_step, const btVector3& pos)
{
    btVector3 currBaseVel;

    if ( m_step == 1 )
        m_prevBaseVel = (pos - m_prevBasePos) / deltaTime;

    if ( m_step > 1 ) {
        currBaseVel = (pos - m_prevBasePos) / deltaTime;
        btVector3 currBaseAcc = (currBaseVel - m_prevBaseVel) / deltaTime;
        for (int i = 0; i < m_numLinks; ++i) {
            btVector3 orient = m_multiBody->getLink(i).m_dVector;
            orient = m_multiBody->localDirToWorld(i, orient);
            btVector3 torque = btCross(orient, -currBaseAcc * m_multiBody->getLink(i).m_mass / (i + 1));
            // TODO at the beginning, the force is too big. but at the next Cycle, it is small. NEED to check and fix.
            torque *= m_linearDragEffect;
            if (torque.length() > MAX_DRAG_FORCE) {
                printf("drag force %f is too large\n", torque.length());
                torque = torque.normalized() * MAX_DRAG_FORCE;
            }
            // TODO addlink is different from addJointTorqueMultiDof, and the later looks better in animation.
//        m_multiBody->addLinkTorque(i, torque);
            for (int d = 0; d < m_multiBody->getLink(i).m_dofCount; d++) {
                m_multiBody->addJointTorqueMultiDof(i, d, torque[d]);
            }
        }
        m_prevBaseVel = currBaseVel;
    }
}

void Skeleton::applyBaseCentrifugalForce(float deltaTime, int m_step, const btVector3& pos)
{

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
//    applySpringForce(deltaTime);
//
////    // calculate the drag force
//    btTransform trans = transformBase(m_time, deltaTime);  // in the beginning, m_time == 0

//    applyBaseCentrifugalForce(deltaTime, m_step, btTransform());
//    m_prevBasePos = trans;

    btVector3 pos = m_move.getPos(m_time, deltaTime);
//    applyBaseLinearDragForce(deltaTime, m_step, pos);
    m_prevBasePos = pos;

    // p2p
    m_p2p->setPivotInB(pos);

    moveCollider(pos);

    // capture the frames
//    {
//        const char* gPngFileName = "multibody";
//        sprintf(mfileName, "%s_%d.png", gPngFileName, m_step);
//        this->m_guiHelper->getAppInterface()->dumpNextFrameToPng(mfileName);
//    }

    // step and update positions
    btScalar  fixedTimeStep = deltaTime / btScalar(2.0);
    m_dynamicsWorld->stepSimulation(deltaTime, 2, fixedTimeStep);
    if ( WIRE_FRAME ) {
        m_dynamicsWorld->debugDrawWorld();
    }

    m_time += deltaTime;
    m_step += 1;
}


class CommonExampleInterface* Dof6Spring2CreateFunc(CommonExampleOptions& options)
{
//    return new Dof6Spring2Setup(options.m_guiHelper);
	return new Skeleton(options.m_guiHelper);
}
