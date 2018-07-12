#include "rtree.h"
#include <rtreenode.h>
#include <utility>
#include <list>
#include <iostream>
#include <chrono>

double rTree::weightedNorm(Eigen::VectorXd x, int dimension){
    double norm(0);

    for(int i = 0; i < dimension; i++){
        norm += weights[i]*x[i]*x[i];
    }

    return std::sqrt(norm);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> rTree::enlargeMBR(rTreeNode *entry1 , rTreeNode *entry2){
    Eigen::VectorXd Snew(entry1->S), Tnew(entry1->T);

    for(int i = 0 ; i < 3; i++){
        if(entry2->S[i] < entry1->S[i])
            Snew[i] = entry2->S[i];
        if(entry2->T[i] > entry1->T[i])
            Tnew[i] = entry2->T[i];
    }

    return(std::make_pair(Snew, Tnew));
}

double rTree::minEnlargement(rTreeNode *entry1 , rTreeNode *entry2){
    std::pair<Eigen::VectorXd, Eigen::VectorXd> newDVertices(enlargeMBR(entry1, entry2));

    return weights[0]*weights[1]*weights[2]*(std::abs(newDVertices.first[0] - newDVertices.second[0])*std::abs(newDVertices.first[1] - newDVertices.second[1])*std::abs(newDVertices.first[2] - newDVertices.second[2]) - std::abs(entry1->S[0] - entry1->T[0])*std::abs(entry1->S[1] - entry1->T[1])*std::abs(entry1->S[2] - entry1->T[2]));
}

rTreeNode *rTree::chooseSubtree(rTreeNode *rootNode, rTreeNode *entry){
    if((rootNode->nodeType)){
        return rootNode;
    }

    rTreeNode *chosenNode(*((rootNode->children).begin()));
    std::list<rTreeNode *>::iterator nodesIterator((rootNode->children).begin());
    double minExpansion(minEnlargement(*nodesIterator, entry));

    nodesIterator++;
    while(nodesIterator != (rootNode->children).end()){
        double newMinExpansion(minEnlargement(*nodesIterator, entry));

        if(newMinExpansion < minExpansion){
            minExpansion = newMinExpansion;
            chosenNode = *nodesIterator;
        }
        else if(newMinExpansion == minExpansion){
            if(std::abs(((*nodesIterator)->S)[0] - ((*nodesIterator)->T)[0])*std::abs(((*nodesIterator)->S)[1] - ((*nodesIterator)->T)[1])*std::abs(((*nodesIterator)->S)[2] - ((*nodesIterator)->T)[2]) <= std::abs((chosenNode->S)[0] - (chosenNode->T)[0])*std::abs((chosenNode->S)[1] - (chosenNode->T)[1])*std::abs((chosenNode->S)[2] - (chosenNode->T)[2]))
                chosenNode = *nodesIterator;
        }
        nodesIterator++;
    }

    return chooseSubtree(chosenNode, entry);
}

std::pair<rTreeNode*, rTreeNode*> rTree::splitNode(rTreeNode *node, rTreeNode *entry){
        double d(0);
        rTreeNode *newNode;
        std::list<rTreeNode *> entries(node->children);
        std::list<rTreeNode *>::iterator it1(entries.begin()), it2(entries.begin());
        std::list<rTreeNode *>::iterator seed1(entries.begin()), seed2(entries.begin());
        std::pair<Eigen::VectorXd, Eigen::VectorXd> dVertices;

        //PICK SEEDS
        node->children.clear();
        entries.emplace_back(entry);
        while(it1 != entries.end()){
            while(it2 != entries.end()){
                if(it1 != it2){
                    std::pair<Eigen::VectorXd, Eigen::VectorXd> tVertices(enlargeMBR(*it1, *it2));

                    double dNew(weights[0]*weights[1]*weights[2]*(std::abs(tVertices.first[0] - tVertices.second[0])*std::abs(tVertices.first[1] - tVertices.second[1])*std::abs(tVertices.first[2] - tVertices.second[2]) - std::abs((*it1)->S[0] - (*it1)->T[0])*std::abs((*it1)->S[1] - (*it1)->T[1])*std::abs((*it1)->S[2] - (*it1)->T[2]) - std::abs((*it2)->S[0] - (*it2)->T[0])*std::abs((*it2)->S[1] - (*it2)->T[1])*std::abs((*it2)->S[2] - (*it2)->T[2])));
                    if(d < dNew){
                        d = dNew;
                        seed1 = it1;
                        seed2 = it2;
                    }
                }
                it2++;
            }
            it1++;
        }

        node->children.emplace_back(*seed1);
        (*seed1)->parentNode = node;
        node->S = (*seed1)->S;
        node->T = (*seed1)->T;
        newNode = new rTreeNode(node->nodeType, m, M, Eigen::VectorXd((*seed2)->S), Eigen::VectorXd((*seed2)->T), nullptr);
        newNode->children.emplace_back(*seed2);
        (*seed2)->parentNode = newNode;
        entries.erase(seed1);
        entries.erase(seed2);

        //CHECK
        while(entries.size() != 0){
            rTreeNode *selectedNode = newNode;

            it1 = entries.begin();
            it2 = entries.begin();

            if((entries.size() + node->children.size() == m) || (entries.size() + newNode->children.size() == m)){
                if(entries.size() + node->children.size() == m)
                    selectedNode = node;

                while(entries.size() != 0){
                    it1 = entries.begin();
                    (*it1)->parentNode = selectedNode;
                    selectedNode->children.emplace_back(*it1);
                    dVertices = enlargeMBR(selectedNode, *it1);
                    selectedNode->S = dVertices.first;
                    selectedNode->T = dVertices.second;
                    entries.erase(it1);
                }

                break;
            }

            //PICK NEXT
            double maxDifference(0);
            double difference;

            while(it1 != entries.end()){
                difference = std::abs(minEnlargement(node, *it1) - minEnlargement(newNode, *it1));
                if (maxDifference < difference){
                    maxDifference = difference;
                    it2 = it1;
                }
                it1++;
            }
            //nasli najbolju

            double enlargement(minEnlargement(node, *it2) - minEnlargement(newNode, *it2));

            if(enlargement < 0){
                selectedNode = node;
            }
            else if(enlargement == 0){
                double areaDifference(weights[0]*weights[1]*weights[2]*(std::abs(node->S[0] - node->T[0])*std::abs(node->S[1]- node->T[1])*std::abs(node->S[2]- node->T[2]) - std::abs(newNode->S[0] - newNode->T[0])*std::abs(newNode->S[1] - newNode->T[1])*std::abs(newNode->S[2] - newNode->T[2])));

                if(areaDifference < 0){
                    selectedNode = node;
                }
                else if(areaDifference == 0){
                    if(node->children.size() <= newNode->children.size()){
                        selectedNode = node;
                    }
                }
            }

            (*it2)->parentNode = selectedNode;
            selectedNode->children.emplace_back(*it2);
            dVertices = enlargeMBR(selectedNode, *it2);
            selectedNode->S = dVertices.first;
            selectedNode->T = dVertices.second;
            entries.erase(it2);
        }

        return std::make_pair(node, newNode);
}

void rTree::adjustTree(rTreeNode *modifiedNode, rTreeNode *newNode){
    if(modifiedNode->parentNode == nullptr){
        if(newNode != nullptr){
            depth++;
            rTreeRoot = new rTreeNode(0, m, M, maxS, maxT, nullptr);
            rTreeRoot->children.emplace_back(modifiedNode);
            rTreeRoot->children.emplace_back(newNode);
            modifiedNode->parentNode = rTreeRoot;
            newNode->parentNode = rTreeRoot;
        }
        return;
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd> dVertices(enlargeMBR(modifiedNode->parentNode, modifiedNode));

    modifiedNode->parentNode->S = dVertices.first;
    modifiedNode->parentNode->T = dVertices.second;

    if(newNode != nullptr){
        if(modifiedNode->parentNode->children.size() < M){
            modifiedNode->parentNode->children.emplace_back(newNode);
            newNode->parentNode = modifiedNode->parentNode;
            dVertices = enlargeMBR(modifiedNode->parentNode, newNode);
            modifiedNode->parentNode->S = dVertices.first;
            modifiedNode->parentNode->T = dVertices.second;
            adjustTree(modifiedNode->parentNode, nullptr);
            return;
        }

        std::pair<rTreeNode*, rTreeNode*> splitPair(splitNode(modifiedNode->parentNode, newNode));

        adjustTree(splitPair.first, splitPair.second);
        return;
    }

    adjustTree(modifiedNode->parentNode, nullptr);
}

rTreeNode* rTree::insertPoint(rTreeNode *rootNode, Eigen::VectorXd S, Eigen::VectorXd T, rrtNode *treeNode){
    rTreeNode *node(new rTreeNode(1, m, M, S, T, treeNode));
    rTreeNode *insertionPoint = chooseSubtree(rootNode, node);
    size++;

    if(insertionPoint->children.size() < M){
        std::pair<Eigen::VectorXd, Eigen::VectorXd> dVertices(enlargeMBR(insertionPoint, node));

        node->parentNode = insertionPoint;
        insertionPoint->children.emplace_back(node);
        insertionPoint->S = dVertices.first;
        insertionPoint->T = dVertices.second;
        adjustTree(insertionPoint, nullptr);
    }
    else{
        std::pair<rTreeNode*, rTreeNode*> splitPair(splitNode(insertionPoint, node));
        adjustTree(splitPair.first, splitPair.second);
    }

    return node;
}

double rTree::minDist(rTreeNode *point, rTreeNode *boundingBox){
    Eigen::VectorXd R(point->S);

    for(int i = 0; i < 2; i++){
        if(point->S[i] < boundingBox->S[i])
            R[i] = boundingBox->S[i];
        else if(point->S[i] > boundingBox->T[i])
            R[i] = boundingBox->T[i];
    }

    return weightedNorm(point->S - R, weights.size());
}

void rTree::kNearestNeighborSearch(double k, rTreeNode *node, rTreeNode *qPoint, std::list<rTreeNode *> &nearestNeighbors, double &nearestDistance){
    std::list<rTreeNode *>::iterator childrenIterator(node->children.begin());
    std::list<rTreeNode *> activeBranchList;//{*childrenIterator};
    std::list<rTreeNode *>::iterator ablIterator, ablIterator2;

    if(k == 0)
            return;

    if(node->nodeType){
        while(childrenIterator != node->children.end()){
            double distance(weightedNorm((*childrenIterator)->S - qPoint->S, weights.size()));
            if((distance < k)){
                nearestNeighbors.emplace_front(*childrenIterator);
            }
            childrenIterator++;
        }
        return;
    }

    activeBranchList = node->children;
    activeBranchList.sort([this, qPoint](rTreeNode *n1,rTreeNode *n2){return minDist(qPoint, n1) <= minDist(qPoint, n2);});
    ablIterator = activeBranchList.begin();

    while(ablIterator != activeBranchList.end()){
        bool innerBreak(false), ablBool(false);
        ablIterator2 = ablIterator;

        while(ablIterator2 != activeBranchList.end()){
            if((minDist(qPoint,*ablIterator2) >= k)){
                ablBool = (ablIterator == ablIterator2);
                activeBranchList.erase(ablIterator2);
                if(ablBool)
                    ablIterator--;

                innerBreak = true;
                break;
            }
            ablIterator2++;
        }
        if(innerBreak)
            continue;

        kNearestNeighborSearch(k, *ablIterator, qPoint, nearestNeighbors, nearestDistance);
        ablIterator++;
    }
}

void rTree::nearestNeighborSearch(rTreeNode *node, rTreeNode *qPoint, rTreeNode *&nearestNeighbor, double &nearestDistance){
    std::list<rTreeNode *>::iterator childrenIterator(node->children.begin());
    std::list<rTreeNode *> activeBranchList{*childrenIterator};
    std::list<rTreeNode *>::iterator ablIterator, ablIterator2;

    if(node->nodeType){
        while(childrenIterator != node->children.end()){
            double distance(weightedNorm((*childrenIterator)->S - qPoint->S, weights.size()));
            if((distance < nearestDistance)){
                nearestDistance = distance;
                nearestNeighbor = (*childrenIterator);
            }
            childrenIterator++;
        }
        return;
    }

    activeBranchList = node->children;
    activeBranchList.sort([this, qPoint](rTreeNode *n1 ,rTreeNode *n2){return minDist(qPoint, n1) <= minDist(qPoint, n2);});
    ablIterator = activeBranchList.begin();

    while(ablIterator != activeBranchList.end()){
        bool innerBreak(false), ablBool(false);
        ablIterator2 = ablIterator;

        while(ablIterator2 != activeBranchList.end()){
            if((minDist(qPoint, *ablIterator2) >= nearestDistance)){
                ablBool = (ablIterator == ablIterator2);
                activeBranchList.erase(ablIterator2);

                if(ablBool)
                    ablIterator--;

                innerBreak = true;
                break;
            }
            ablIterator2++;
        }
        if(innerBreak)
            continue;

        nearestNeighborSearch(*ablIterator, qPoint, nearestNeighbor, nearestDistance);
        ablIterator++;
    }
}

rTree::rTree(int _m, int _M, Eigen::VectorXd dim1, Eigen::VectorXd dim2, std::vector<double> wghts): m(_m), M(_M), maxS(dim1), maxT(dim2), depth(0), size(0), rTreeRoot(new rTreeNode(1, _m, _M, dim1, dim2, nullptr)), weights(wghts){}
