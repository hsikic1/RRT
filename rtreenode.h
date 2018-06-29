#ifndef RTREENODE_H
#define RTREENODE_H
#include <Eigen/Dense>
#include <list>
#include <vector>
#include <rrtnode.h>

class rTreeNode
{   
public:
    unsigned int m,M;
    bool nodeType;
    bool root;
    rTreeNode *parentNode;
    Eigen::VectorXd S, T;
    std::list<rTreeNode *> children;
    rrtNode *treeNode;
    //std::list<Eigen::VectorXd> points;

    rTreeNode(bool type, int _m, int _M, Eigen::VectorXd s, Eigen::VectorXd t, rrtNode *node);
    ~rTreeNode();
};

#endif // RTREENODE_H
