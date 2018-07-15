#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <rtreenode.h>
#include <rtree.h>
#include <rrtnode.h>
#include <vector>
#include <iostream>
#include <string>
#include <QMainWindow>
//#include <gtest/gtest.h>
#include "fcl/math/bv/utility.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/detail/gjk_solver_indep.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    int sum;
    int segNumber;
    std::vector<double> weights;
    QVector<double> x, y;
    std::vector<fcl::Vector3<double>> environmentVertices;
    std::vector<fcl::Triangle> environmentTriangles;
    std::vector<std::vector<fcl::Vector3<double>>> robotSegVertices;
    std::vector<std::vector<fcl::Triangle>> robotSegTriangles;
    std::vector<Eigen::VectorXd> fKine;
    std::vector<double> segLengths;
    rrtNode *lastNode;

    void traverseTree(rTreeNode *node);
    bool simpleCollisionCheck(Eigen::Vector2d center, double radius, Eigen::VectorXd A, Eigen::VectorXd B);
    double weightedNorm(Eigen::VectorXd x, int dimension);
    void makePlot();

private slots:
    void parseSTL(std::string path, std::vector<fcl::Vector3<double>> &vertices, std::vector<fcl::Triangle> &triangles);
    void loadMeshes(std::string robotPath, std::string envPath);
    void initRRT(Eigen::VectorXd xinit, Eigen::VectorXd xgoal, std::list<std::pair<Eigen::Vector2d , double>> &obstacles);
    rrtNode *extendRRTStar(rrtNode *&epsilonNode, rTree *nodes, Eigen::VectorXd xrandom, std::list<std::pair<Eigen::Vector2d, double>> &obstacles, Eigen::VectorXd xgoal);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
