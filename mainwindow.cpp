#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <Eigen/Dense>
#include <rtree.h>
#include <qcustomplot.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
//#include <gtest/gtest.h>
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl/common/unused.h"
#include "fcl/math/constants.h"
#include "fcl/math/triangle.h"

#define epsilon 0.2
#define R 10.5959179423
#define segs 3
#define a1 2.2
#define a2 1.7
#define a3 0.9


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->customPlot->xAxis->setRange(-5,2*M_PI+0.1);
    ui->customPlot->yAxis->setRange(-5,2*M_PI+0.1);

    weights = {1.0000,    0.6833,    0.3377};
    segNumber = 3;

    Eigen::Vector2d xgoal(3.9,2.2);
    Eigen::VectorXd xinit(segNumber);
    xinit(0) = M_PI;
    xinit(1) = 0;
    xinit(2) = 0;
    segLengths = {2.2, 1.7, 0.9};

    std::list<std::pair<Eigen::Vector2d, double>> obstacles{std::make_pair(Eigen::Vector2d(-1.7*std::sqrt(2),1.7*sqrt(2)), 1),std::make_pair(Eigen::Vector2d(1.7*std::sqrt(2),1.7*sqrt(2)), 0.9),std::make_pair(Eigen::Vector2d(1.7*std::sqrt(2),-0.5), 0.8)};
    initRRT(xinit, xgoal, obstacles);

    //loadMeshes("/home/haris/Desktop/robot/solid.robot", "");
}

void MainWindow::parseSTL(std::string path, std::vector<fcl::Vector3<double>> &vertices, std::vector<fcl::Triangle> &triangles){
    std::ifstream fileStream(path, std::ios::binary);
    unsigned char dummyByte(0);
    unsigned char countUint8[4];
    unsigned char atributeUint8[2];
    unsigned int meshCount(0);

    if(!fileStream)
        std::cout << "Nema fajla";

    for(int i = 0; i < 80; i++){
        fileStream.read(reinterpret_cast<char *>(&dummyByte), sizeof dummyByte);
    }

    fileStream.read(reinterpret_cast<char *>(countUint8), sizeof countUint8);
    for(int i = 0; i < 4; i++){
        meshCount |= countUint8[i] << i;
    }

    for(unsigned int i = 0; i < meshCount; i++){
        for(int i = 0; i < 4; i++){
            float x,y,z;
            fcl::Triangle tempTriangle(vertices.size(), vertices.size() + 1, vertices.size() + 2);

            fileStream.read(reinterpret_cast<char *>(&x), sizeof x);
            fileStream.read(reinterpret_cast<char *>(&y), sizeof y);
            fileStream.read(reinterpret_cast<char *>(&z), sizeof z);
            vertices.emplace_back(fcl::Vector3<double>(x,y,z));
            triangles.emplace_back(tempTriangle);
        }
        fileStream.read(reinterpret_cast<char *>(&atributeUint8), sizeof atributeUint8);
    }
}

void MainWindow::loadMeshes(std::string robotPath, std::string envPath){
    std::ifstream jsonFileStream(robotPath);
    rapidjson::IStreamWrapper jsonStreamWrapper(jsonFileStream);
    rapidjson::Document docInst;

    docInst.ParseStream(jsonStreamWrapper);
    parseSTL("/home/haris/Desktop/robot/environment.stl", environmentVertices, environmentTriangles);

    for(unsigned int k = 0; k < docInst.MemberCount(); k++){
        std::ostringstream tempStream("seg", std::ios::ate);
        tempStream << k + 1;

        robotSegVertices.push_back(std::vector<fcl::Vector3<double>>());
        robotSegTriangles.push_back(std::vector<fcl::Triangle>());

        parseSTL(std::string("/home/haris/Desktop/robot/") + docInst[tempStream.str().c_str()].GetObject()["parts"].GetString(), robotSegVertices[k], robotSegTriangles[k]);
    }
}

void MainWindow::traverseTree(rTreeNode *node){
    std::list<rTreeNode *>::iterator nIter(node->children.begin());

    auto line1 = new QCPItemLine(ui->customPlot);
    auto line2 = new QCPItemLine(ui->customPlot);
    auto line3 = new QCPItemLine(ui->customPlot);
    auto line4 = new QCPItemLine(ui->customPlot);

    line1->start->setCoords(node->T[0], node->T[1]);
    line2->start->setCoords(node->T[0], node->T[1]);
    line1->end->setCoords(node->T[0], node->S[1]);
    line2->end->setCoords(node->S[0], node->T[1]);

    line3->start->setCoords(node->S[0], node->S[1]);
    line4->start->setCoords(node->S[0], node->S[1]);
    line3->end->setCoords(node->S[0], node->T[1]);
    line4->end->setCoords(node->T[0], node->S[1]);

    if(node->children.back()->children.size() == 0){
        auto piter(node->children.begin());
        while(piter != node->children.end()){
            if((*piter)->treeNode->parentNode != nullptr){
                auto line = new QCPItemLine(ui->customPlot);
                line->start->setCoords((*piter)->treeNode->nodePosition[0], (*piter)->treeNode->nodePosition[1]);
                line->end->setCoords((*piter)->treeNode->parentNode->nodePosition[0], (*piter)->treeNode->parentNode->nodePosition[1]);
                line->setHead(QCPLineEnding(QCPLineEnding::esDisc, 3));
                line->setTail(QCPLineEnding(QCPLineEnding::esDisc, 3));

                QPen pen;
                pen.setWidth(1);
                pen.setColor(QColor(66,131,244));
                line->setPen(pen);
            }
            piter++;
        }
        sum += node->children.size();
        return;
    }

    while(nIter != node->children.end()){
        traverseTree(*nIter);
        nIter++;
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

double MainWindow::weightedNorm(Eigen::VectorXd x, int dimension){
    double norm(0);

    for(int i = 0; i < dimension; i++){
        norm += weights[i]*x[i]*x[i];
    }

    return std::sqrt(norm);
}

bool MainWindow::simpleCollisionCheck(Eigen::Vector2d center, double radius, Eigen::VectorXd A, Eigen::VectorXd B){
    //QVector<double>x, y;
    fKine = std::vector<Eigen::VectorXd>(segNumber + 1, Eigen::VectorXd(2));
    fKine[0](0) = 0;
    fKine[0](1) = 0;
    /*Eigen::VectorXd fk1(2);
    Eigen::VectorXd fk2(2);
    Eigen::VectorXd fk3(2);*/

    for(int i = 10; i > 0; i--){
        Eigen::VectorXd step(A + (i/10.0)*epsilon*(B-A).normalized());
        std::vector<double> angles(segNumber, 0);
        //double stepSum(0);

        for(int j = 1; j < segNumber + 1; j++){
            fKine[j](0) = 0;
            fKine[j](1) = 0;
            angles[j-1] += step[j - 1];
            if(j - 1 > 0)
                angles[j-1] += angles[j-2];

            for(int k = 0; k < j; k++){
                fKine[j](0) += segLengths[k]*cos(angles[k]);
                fKine[j](1) += segLengths[k]*sin(angles[k]);
            }

            if((((center - (fKine[j] - fKine[j-1] - ((center-fKine[j-1]).dot((fKine[j]-fKine[j-1]).normalized()))*((fKine[j]-fKine[j-1]).normalized()))).norm() <= radius) || ((center-fKine[j-1]).norm() <= radius) || ((center-fKine[j]).norm() <= radius)))
                return true;
        }
    }

    return false;
}

//8 16
void MainWindow::initRRT(Eigen::VectorXd xinit, Eigen::Vector2d xgoal, std::list<std::pair<Eigen::Vector2d , double>> &obstacles){
    Eigen::VectorXd S(segNumber), T(segNumber);
    for(int i = 0; i < segNumber; i++){
        S(i) = 0;
        T(i) = 2*M_PI;
    }

    QPen pen;
    rTree *nodes(new rTree(8,16,S, T, weights));
    rrtNode *initRRTN(new rrtNode(xinit, nullptr, 0));
    rrtNode *rrtPtr(nullptr), *rrtNew(nullptr);

    nodes->insertPoint(nodes->rTreeRoot, xinit, xinit, initRRTN);
    nodes->initNode = new rTreeNode(1,2,5,xinit,xinit, initRRTN);

    auto t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < 25000; i++){
        Eigen::VectorXd xrandom(segNumber);

        for(int j = 0; j < segNumber; j++){
            std::random_device rd1{};
            std::seed_seq seed1{rd1(), rd1(), rd1(), rd1(), rd1(), rd1(), rd1(), rd1()};
            std::mt19937 engine1(seed1);
            std::uniform_real_distribution<double> dist1{0.0, 2*M_PI};
            xrandom(j) = dist1(engine1);
        }

        rrtNew = extendRRTStar(nodes, xrandom, obstacles, xgoal);
        if(rrtNew != nullptr)
            rrtPtr = rrtNew;
        if(rrtPtr != nullptr)
            break;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << duration << std::endl;

    //traverseTree(nodes->rTreeRoot);
    /*if(rrtPtr == nullptr)
        return;

    while(rrtPtr->parentNode != nullptr){
        auto line = new QCPItemLine(ui->customPlot);
        line->start->setCoords(rrtPtr->nodePosition[0], rrtPtr->nodePosition[1]);
        line->end->setCoords(rrtPtr->parentNode->nodePosition[0], rrtPtr->parentNode->nodePosition[1]);
        line->setHead(QCPLineEnding(QCPLineEnding::esDisc, 3));
        line->setTail(QCPLineEnding(QCPLineEnding::esDisc, 3));
        QPen pen;
        pen.setWidth(2);
        pen.setColor(QColor(163, 23, 13));
        line->setPen(pen);
        rrtPtr = rrtPtr->parentNode;
    }*/

    if(rrtPtr == nullptr){
        return;
    }

    ui->customPlot->addGraph();
    while(rrtPtr != nullptr){
        std::vector<double> angles(segNumber, 0);
        Eigen::VectorXd fKineVec1(2);
        Eigen::VectorXd fKineVec2(2);

        for(int i = 1; i < segNumber + 1; i++){
            fKineVec1(0) = 0;
            fKineVec1(1) = 0;
            fKineVec2(0) = 0;
            fKineVec2(1) = 0;
            angles[i -1] += rrtPtr->nodePosition[i -1];

            if(i - 1 > 0)
                angles[i - 1] += angles[i -2];

            for(int j = 0; j < i; j++){
                fKineVec2(0) += segLengths[j]*std::cos(angles[j]);
                fKineVec2(1) += segLengths[j]*std::sin(angles[j]);
                if(j < i - 1){
                    fKineVec1(0) += segLengths[j]*std::cos(angles[j]);
                    fKineVec1(1) += segLengths[j]*std::sin(angles[j]);
                }
            }

            auto line = new QCPItemLine(ui->customPlot);

            if(rrtPtr == rrtNew){
                pen.setWidth(5);
                pen.setColor(QColor(66,131,244));
                line->setPen(pen);
            }

            if(rrtPtr->parentNode == nullptr){
                pen.setWidth(5);
                pen.setColor(QColor(163, 23, 13));
                line->setPen(pen);

            }

            line->start->setCoords(fKineVec1[0], fKineVec1[1]);
            line->end->setCoords(fKineVec2[0], fKineVec2[1]);
        }

        /*line->start->setCoords(rrtPtr->nodePosition[0], rrtPtr->nodePosition[1]);
        line->end->setCoords(rrtPtr->parent->nodePosition[0], rrtPtr->parent->nodePosition[1]);
        line->setHead(QCPLineEnding(QCPLineEnding::esDisc, 3));
        line->setTail(QCPLineEnding(QCPLineEnding::esDisc, 3));
        QPen pen;
        pen.setWidth(2);
        pen.setColor(QColor(163, 23, 13));
        line->setPen(pen);*/

        /*if(rrtPtr->parentNode == nullptr){
            pen.setWidth(5);
            pen.setColor(QColor(163, 23, 13));
            line->setPen(pen);
            line1->setPen(pen);
            line2->setPen(pen);
            line3->setPen(pen);
            break;
        }*/
        rrtPtr = rrtPtr->parentNode;
    }

    std::list<std::pair<Eigen::Vector2d, double>>::iterator obstacleIterator(obstacles.begin());

    QVector<double> x(0), y(0);
    ui->customPlot->addGraph();

    while(obstacleIterator != obstacles.end()){
        for(int i = 0; i < 100; i++){
            x.append((*obstacleIterator).second*std::cos(i*M_PI/50) + (*obstacleIterator).first[0]);
            y.append((*obstacleIterator).second*std::sin(i*M_PI/50) + (*obstacleIterator).first[1]);
        }
        ui->customPlot->graph()->addData(x,y);
        x.clear();
        y.clear();
        obstacleIterator++;
        ui->customPlot->addGraph();
    }
}

//EXTEND RRT*
rrtNode *MainWindow::extendRRTStar(rTree *nodes, Eigen::VectorXd xrandom, std::list<std::pair<Eigen::Vector2d , double>> &obstacles, Eigen::Vector2d xgoal){
    std::list<std::pair<Eigen::Vector2d, double>>::iterator obstacleIterator(obstacles.begin());
    std::list<rTreeNode *>::iterator nIterator;
    std::list<rTreeNode *> nearestNeighbors;
    rTreeNode *xmin(nodes->initNode);
    double dist(weightedNorm(xrandom - nodes->initNode->S, segNumber)); //norm()

    nodes->nearestNeighborSearch(nodes->rTreeRoot, new rTreeNode(1,2,5, xrandom, xrandom, nullptr), xmin, dist); //

    Eigen::VectorXd xepsilon(epsilon*(xrandom - xmin->S).normalized() + xmin->S); //srediti normiranje

    while(obstacleIterator != obstacles.end()){
        if(simpleCollisionCheck(obstacleIterator->first, obstacleIterator->second, xmin->S, xepsilon))//
            return nullptr;
        obstacleIterator++;
    }

    double m = std::min(2*epsilon, R*std::sqrt(std::log10(nodes->size)/((double)(nodes->size))));
    dist = weightedNorm(xrandom - nodes->initNode->S,segNumber);  //.norm();
    nodes->kNearestNeighborSearch(m, nodes->rTreeRoot, new rTreeNode(1,2,5, xepsilon, xepsilon, nullptr), nearestNeighbors, dist);

    nIterator = nearestNeighbors.begin();
    while(nIterator != nearestNeighbors.end()){
        bool breakFlag(false);

        obstacleIterator = obstacles.begin();
        while(obstacleIterator != obstacles.end()){
            if(simpleCollisionCheck(obstacleIterator->first, obstacleIterator->second, (*nIterator)->S, xepsilon)){
                breakFlag = true;
                break;
            }
            obstacleIterator++;
        }

        if(breakFlag){
            nearestNeighbors.erase(nIterator);
            nIterator--;
            continue;
        }

        if((*nIterator)->treeNode->cost + weightedNorm((*nIterator)->S - xepsilon, segNumber) < xmin->treeNode->cost + weightedNorm(xepsilon - xmin->S, segNumber)){
            xmin = *nIterator;
        }
        nIterator++;
    }

    rrtNode *insertedNode(new rrtNode(xepsilon, xmin->treeNode, xmin->treeNode->cost + weightedNorm(xepsilon - xmin->S, segNumber)));
    nodes->insertPoint(nodes->rTreeRoot, xepsilon, xepsilon, insertedNode);

    nIterator = nearestNeighbors.begin();
    while(nIterator != nearestNeighbors.end()){
        if((*nIterator != xmin) && (*nIterator)->treeNode->cost > weightedNorm((*nIterator)->S - xepsilon, segNumber) + insertedNode->cost){
            (*nIterator)->treeNode->parentNode = insertedNode;
            (*nIterator)->treeNode->cost = weightedNorm((*nIterator)->S - xepsilon, segNumber) + insertedNode->cost;
        }
        nIterator++;
    }

    double confSum(0);
    Eigen::VectorXd fk(2);
    fk(0) = 0;
    fk(1) = 0;

    for(int i = 0; i < segNumber; i++){
        confSum += xepsilon[i];
        fk(0) += segLengths[i]*std::cos(confSum);
        fk(1) += segLengths[i]*std::sin(confSum);
    }

    if((fk - xgoal).norm() <= 0.1)
        return insertedNode;

    return nullptr;
}
