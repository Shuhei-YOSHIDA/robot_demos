#pragma once
#include <Eigen/Core>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/Jacobian.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/EulerIntegration.h>
#include "pseudoInverse.h"
using namespace std;
using namespace Eigen;

class Task {// interface class

public:
    virtual ~Task(){};
    virtual VectorXd g(rbd::MultiBody mb, rbd::MultiBodyConfig mbc) = 0;
    virtual MatrixXd J(rbd::MultiBody mb, rbd::MultiBodyConfig mbc) = 0;
    string _name;
    string _type;
};

typedef boost::shared_ptr<Task> TaskPtr;
typedef std::vector<pair<double, TaskPtr>> MultiTaskPtr;//priority and Task

class BodyTask : public Task {

public:
    BodyTask(rbd::MultiBody mb, string bodyName, sva::PTransformd X_O_T, 
             sva::PTransformd X_b_p = sva::PTransformd::Identity(), string part = "",  string name = "BodyTask")
    {
        /*
        Compute the error and the jacobian to target a static frame for a body.
        Parameters: 
        - mb: MultiBody
        - bodyId: ID of the body that should move
        - X_0_T: targeted frame (PTransformd)
        - X_b_p: static frame on the body bodyId
        - part: "position", "rotation", something else
        */
        _name = name;
        _bodyName = bodyName;
        _bodyIndex = mb.bodyIndexByName(bodyName);
        _X_O_T = X_O_T;
        _X_b_p = X_b_p;
        _jac = rbd::Jacobian(mb, bodyName);
        _jac_mat_sparse = MatrixXd::Zero(6, mb.nrDof());
        _type = "BodyTask";
        _part = part;
    }

    sva::PTransformd X_O_p(rbd::MultiBodyConfig mbc) 
    {
        sva::PTransformd X_O_b(mbc.bodyPosW[_bodyIndex]);
        sva::PTransformd xop = _X_b_p * X_O_b;
        return xop; 
    }

    virtual VectorXd g(rbd::MultiBody mb, rbd::MultiBodyConfig mbc)
    {
        auto X_O_p = this->X_O_p(mbc);
        auto g_body = sva::transformError(_X_O_T, X_O_p);//MotionVec
        if      (_part == "position") return g_body.vector().segment(3, 3);
        else if (_part == "rotation") return g_body.vector().segment(0, 3);
        return g_body.vector();
    }

    virtual MatrixXd J(rbd::MultiBody mb, rbd::MultiBodyConfig mbc)
    {
        sva::PTransformd X_O_p = this->X_O_p(mbc);
        // Set transformation in Origin orientation frame
        sva::PTransformd X_O_p_0 = sva::PTransformd(X_O_p.rotation()).inv() * X_O_p;
        MatrixXd jac_mat_dense = _jac.jacobian(mb, mbc, X_O_p_0);
        _jac.fullJacobian(mb, jac_mat_dense, _jac_mat_sparse);
        if      (_part == "position") return _jac_mat_sparse.block(3, 0, 3, mb.nrDof());
        else if (_part == "rotation") return _jac_mat_sparse.block(0, 0, 3, mb.nrDof());
        return _jac_mat_sparse;
    }
    string _bodyName;
    int _bodyIndex;
    sva::PTransformd _X_O_T;
    sva::PTransformd _X_b_p;
    rbd::Jacobian  _jac;
    MatrixXd _jac_mat_sparse;
    string _part;
};


class PostureTask : public Task {

public:
    PostureTask(rbd::MultiBody mb, rbd::MultiBodyConfig mbc, string name = "PostureTask")
    {
        _name = name;
        _q_T = mbc; //need only one reference articular position vector

        //use joint type: prism, rev and Sperical
        auto isDefine =[](rbd::Joint j) { 
            if (j.type() == rbd::Joint::Prism || 
                j.type() == rbd::Joint::Rev   ||
                j.type() == rbd::Joint::Spherical) return true;
            else return false;
        };
        
        // take back joint and joint index that are define
        auto alljoints = mb.joints();
        int count = 0;
        for (auto itr = alljoints.begin(); itr != alljoints.end(); ++itr) {
            if (isDefine(*itr)) {
                _jointIndex.push_back(count);
                _joints.push_back(*itr);
            }
            count++;
        }
        int nrDof = 0;
        for (auto itr = _joints.begin(); itr != _joints.end(); ++itr) {
            nrDof += itr->dof();
        }

        // initialize g
        _g_mat = VectorXd::Zero(nrDof);

        // initialize the jacobian
        _J_mat = MatrixXd::Zero(nrDof, mb.nrDof());
        unsigned int posInG = 0;
        count = 0;
        for (auto itr = _joints.begin(); itr != _joints.end(); ++itr) {
            unsigned int posInDof = mb.jointPosInDof(_jointIndex[count]);
            _J_mat.block(posInG, posInDof, itr->dof(), itr->dof()) = 
                MatrixXd::Identity(itr->dof(), itr->dof());
            posInG += itr->dof();
            count++;
        }
        
        _type = "PostureTask";
    }

    virtual VectorXd g(rbd::MultiBody mb, rbd::MultiBodyConfig mbc)
    {
        auto q = mbc.q; 
        auto jointConfig = mbc.jointConfig;
        unsigned int posInG = 0;
        int count = 0;
        for (auto itr = _joints.begin(); itr != _joints.end(); ++itr) {
            int jIndex = _jointIndex[count];
            if (itr->type() == rbd::Joint::Prism || itr->type() == rbd::Joint::Rev) {
                VectorXd tmp(itr->dof()); tmp << q[jIndex][0] - _q_T.q[jIndex][0];
                _g_mat.segment(posInG, itr->dof()) = tmp; 
            }
            else if (itr->type() == rbd::Joint::Spherical) {
                auto orid = Quaterniond(_q_T.q[jIndex][0], _q_T.q[jIndex][1], 
                                        _q_T.q[jIndex][2], _q_T.q[jIndex][3]).inverse().matrix();
                _g_mat.segment(posInG, itr->dof()) = 
                    sva::rotationError(orid, jointConfig[jIndex].rotation());
            }

            posInG+=itr->dof();
            count++;
        }
        return _g_mat;
    }

    virtual MatrixXd J(rbd::MultiBody mb, rbd::MultiBodyConfig mbc)
    {
        return _J_mat;
    }

    std::vector<unsigned int> _jointIndex;
    std::vector<rbd::Joint> _joints;
    rbd::MultiBodyConfig _q_T;
    VectorXd _g_mat;
    MatrixXd _J_mat;
};

void oneTaskMin(rbd::MultiBody mb, rbd::MultiBodyConfig &mbc, TaskPtr task,
                double delta = 1.0, unsigned int maxIter = 100, double prec = 1e-8)
{
    auto q = rbd::paramToVector(mb, mbc.q);
    unsigned int iterate = 0;
    bool minimizer = false;
    while (iterate < maxIter && !minimizer) {
        // compute task data
        auto g = task->g(mb, mbc);
        auto J = task->J(mb, mbc);
        
        // compute alpha
        // J*alpha = -g
        //cout << "alpha calc" << endl;
        VectorXd alpha;
        alpha = -PseudoInverse(J, 1e-9)*g;//least square ?
    
        // integrate and run the forward kinematic
        //cout << "integration calc" << endl;
        mbc.alpha = rbd::vectorToDof(mb, alpha);
        rbd::eulerIntegration(mb, mbc, delta);
        rbd::forwardKinematics(mb, mbc);
    
        // take the new q vector
        //cout << "new q vector calc" << endl;
        q = rbd::paramToVector(mb, mbc.q);
    
        //H-infinite norm?
        //cout << "infinite norm calc" << endl;
        auto alphaInf = alpha.lpNorm<Infinity>();
    
        // yield the current state
        cout << "iterate " << iterate << endl;
        //cout << " --------------------------- " << endl;
        //cout << "iterate " << iterate << endl
        //     << "q " << endl << q.transpose() << endl
        //     << "alpha " << alpha.transpose() << endl 
        //     << "alphainf " << alphaInf << endl;
        //cout << "g " << endl << g << endl;
        //cout << "J " << endl << J << endl;
        //cout << " --------------------------- " << endl;
    
        // check if the current alpha is a minimizer
        if (alphaInf < prec) minimizer = true;
        iterate++;
    }
}

void manyTaskMin(rbd::MultiBody mb, rbd::MultiBodyConfig &mbc, MultiTaskPtr tasks, 
                 double delta = 1.0, unsigned int maxIter = 100, double prec = 1e-8)
{
    auto q = rbd::paramToVector(mb, mbc.q);
    unsigned int iterate = 0;
    bool minimizer = false;
    while (iterate < maxIter && !minimizer) {
        // compute task data such as w*g, w*J
        std::vector<VectorXd> gList;
        std::vector<MatrixXd> JList;
        int count = 0;
        int height = 0;
        for (auto itr = tasks.begin(); itr != tasks.end(); ++itr) {
            VectorXd gi = itr->first * itr->second->g(mb, mbc);
            MatrixXd Ji = itr->first * itr->second->J(mb, mbc);

            //test check: undefined value causes crash
            for (double var : vector<double>(gi.data(), gi.data()+gi.size())) {
                if (std::isnan(var)) {
                    cout << "gi nan " << endl;
                    while(1);
                }
            }
            for (double var : vector<double>(Ji.data(), Ji.data()+Ji.size())) {
                if (std::isnan(var)) {
                    cout << "Ji nan " << endl;
                    while(1);
                }
            }
            gList.push_back(gi);
            JList.push_back(Ji);
            height+=gi.size();//J.rows == g.rows
            count++;
        }
        //concatinate gList, JList
        int width = mb.nrDof(); 
        VectorXd g(height);
        MatrixXd J(height, width);
        int len = 0;
        for (auto itr = gList.begin(); itr != gList.end(); ++itr) {
            g.segment(len, itr->size()) = *itr;
            len+=itr->size();
        }
        len = 0;
        for (auto itr = JList.begin(); itr != JList.end(); ++itr) {
            J.block(len, 0, itr->rows(), width) = *itr;
            len+=itr->rows();
        }

        
        // compute alpha
        // J*alpha = -g
        VectorXd alpha;
        alpha = -PseudoInverse(J)*g;//least square ?

        // integrate and run the forward kinematics
        mbc.alpha = rbd::vectorToDof(mb, alpha);
        rbd::eulerIntegration(mb, mbc, delta);
        rbd::forwardKinematics(mb, mbc);

        // take the new q vector
        q = rbd::paramToVector(mb, mbc.q);
    
        //H-infinite norm?
        auto alphaInf = alpha.lpNorm<Infinity>();

        // yield the current state
        //cout << " iterate "  << iterate << endl;
        //cout << " --------------------------- " << endl;
        //cout << "iterate " << iterate << endl
        //     << "q " << endl << q.transpose() << endl
        //     << "alpha " << alpha.transpose() << endl 
        //     << "alphainf " << alphaInf << endl;
        //cout << "g " << endl << g << endl;
        //cout << "J " << endl << J << endl;
        //cout << " --------------------------- " << endl;

        // check if the current alpha is a minimizer
        if (alphaInf < prec) minimizer = true;
        iterate++;
    }
    if (iterate>=maxIter) ROS_DEBUG("max_itr");
}


