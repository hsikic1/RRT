#include "rtreenode.h"
#include <Eigen/Dense>

rTreeNode::rTreeNode(bool type, int _m, int _M, Eigen::VectorXd s, Eigen::VectorXd t, rrtNode *node): children(std::list<rTreeNode*>()), parentNode(nullptr), nodeType(type), S(s), T(t), m(_m), M(_M), treeNode(node){}
rTreeNode::~rTreeNode(){}
