#include "rrtnode.h"

rrtNode::rrtNode(Eigen::VectorXd pos, rrtNode *par, double cos): nodePosition(pos), parentNode(par), cost(cos){}
