#ifndef RTREE_H
#define RTREE_H
#include <rtreenode.h>
#include <rrtnode.h>
#include <utility>
#include <list>
#include <vector>


class rTree
{
private:
    unsigned int m, M;
    std::vector<double> weights;
    std::pair<Eigen::VectorXd, Eigen::VectorXd> enlargeMBR(rTreeNode *entry1 , rTreeNode *entry2);
    double minEnlargement(rTreeNode *entry1, rTreeNode *entry2);
    rTreeNode *chooseSubtree(rTreeNode *rootNode, rTreeNode *entry);
    std::pair<rTreeNode*, rTreeNode*> splitNode(rTreeNode *node, rTreeNode *entry);
    void adjustTree(rTreeNode *modifiedNode, rTreeNode *newNode);
public:
    Eigen::VectorXd maxS, maxT;
    int depth, size;
    rTreeNode *rTreeRoot, *initNode;
    rTree(int _m, int _M, Eigen::VectorXd dim1, Eigen::VectorXd dim2, std::vector<double> wghts);
    rTreeNode *insertPoint(rTreeNode *rootNode, Eigen::VectorXd S, Eigen::VectorXd T, rrtNode *treeNode);
    void nearestNeighborSearch(rTreeNode *node, rTreeNode *qPoint, rTreeNode *&nearestNeighbor, double &nearestDistance);
    void kNearestNeighborSearch(double k, rTreeNode *node, rTreeNode *qPoint, std::list<rTreeNode *> &nearestNeighbors, double &nearestDistance);
    double minDist(rTreeNode *point, rTreeNode *boundingBox);
    double minMaxDist(rTreeNode *point, rTreeNode *boundingBox);
    double weightedNorm(Eigen::VectorXd x, int dimension);
};

#endif // RTREE_H
