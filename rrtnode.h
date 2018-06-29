#ifndef RRTNODE_H
#define RRTNODE_H
#include <Eigen/Dense>
#include <list>

class rTreeNode;

class rrtNode
{
public:
    rrtNode *parentNode;
    std::list<rrtNode *> children;
    Eigen::VectorXd nodePosition;
    double cost;
    rTreeNode *rNode;

    rrtNode(Eigen::VectorXd pos, rrtNode *par, double cos);
};

#endif // RRTNODE_H
