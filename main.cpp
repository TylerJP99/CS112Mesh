#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <igl/decimate.h>

#include <igl/collapse_edge.h>
#include <igl/edges.h>
#include <igl/min_heap.h>

Eigen::MatrixXd V, OV;
Eigen::MatrixXi F, OF;
Eigen::VectorXi j;

igl::min_heap<std::tuple<double, int , int>> Q;
igl::min_heap<std::tuple<double,Eigen::Vector4d,int>> Vhats;
Eigen::Matrix<int, Eigen::Dynamic, 2> E;

void computeOptimalContract(
  std::vector<std::pair<Eigen::MatrixXd, int>> Qs,
  Eigen::Matrix<int, Eigen::Dynamic, 2>& E)
{
  
  Eigen::Vector4d Vhat, I(0,0,0,1);
  Eigen::Matrix4d Q1, Q2, Qhat, cry;


  for (unsigned int i = 0; i < E.rows(); i++)
  {
    int e1 = E(i, 0);
    int e2 = E(i, 1);

    for (unsigned int l = 0; l < Qs.size(); l++)
    {
      if (Qs[l].second == e1)
      {
        Q1 = Qs[l].first;
      }
      else if (Qs[l].second == e2)
      {
        Q2 = Qs[l].first;
      }
    }

    Qhat = Q1 + Q2;

    cry(0,0) = Qhat(0,0);
    cry(0,1) = Qhat(0,1);
    cry(0,2) = Qhat(0,2);
    cry(0,3) = Qhat(0,3);
    cry(1,0) = Qhat(0,1);
    cry(1,1) = Qhat(1,1);
    cry(1,2) = Qhat(1,2);
    cry(1,3) = Qhat(1,3);
    cry(2,0) = Qhat(0,2);
    cry(2,1) = Qhat(1,2);
    cry(2,2) = Qhat(2,2);
    cry(2,3) = Qhat(2,3);
    cry(3,0) = 0;
    cry(3,1) = 0;
    cry(3,2) = 0;
    cry(3,4) = 1;

    Vhat = cry.inverse()*I;
    float cost = Vhat.transpose() * Qhat * Vhat;

    Q.emplace(cost, e1, e2);
    Vhats.emplace(cost, Vhat, 0);
  }
}

std::vector<std::pair<Eigen::MatrixXd, int>> computeQ(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  std::vector<std::pair<Eigen::MatrixXd, int>> Qs;
  Eigen::Vector3d p1,p2,p3,v,w;
  Eigen::Vector4d n;

  for (unsigned int i = 0; i < V.rows(); ++i)
  {
    Eigen::MatrixXd Kp(4,4), Q(4,4);
    for (unsigned int rows = 0; rows < F.rows(); rows++)
    {
      if (i == F(rows, 0))
      {
        p1 = V.row(i);
        p2 = V.row(F(rows, 1));
        p3 = V.row(F(rows, 2));
        v = p2 - p1;
        w = p3 - p1;

        n(0) = (v(1) * w(2)) - (v(2) * w(1)); 
        n(1) = (v(2) * w(0)) - (v(0) * w(2)); 
        n(2) = (v(0) * w(1)) - (v(1) * w(0)); 
        n(3) = n(0) * p1(0) + n(1) * p1(1) + n(2) * p1(2);

        // NORMALIZE 
        Kp = n * n.transpose();
        Q += Kp;
      }
      else if (i == F(rows, 1))
      {
        p1 = V.row(i);
        p2 = V.row(F(rows, 0));
        p3 = V.row(F(rows, 2));
        v = p2 - p1;
        w = p3 - p1;

        n(0) = (v(1) * w(2)) - (v(2) * w(1)); 
        n(1) = (v(2) * w(0)) - (v(0) * w(2)); 
        n(2) = (v(0) * w(1)) - (v(1) * w(0)); 
        n(3) = n(0) * p1(0) + n(1) * p1(1) + n(2) * p1(2);

        Kp = n * n.transpose();
        Q += Kp;
      }
      else if (i == F(rows, 2))
      {
        p1 = V.row(i);
        p2 = V.row(F(rows, 1));
        p3 = V.row(F(rows, 0));
        v = p2 - p1;
        w = p3 - p1;

        n(0) = (v(1) * w(2)) - (v(2) * w(1)); 
        n(1) = (v(2) * w(0)) - (v(0) * w(2)); 
        n(2) = (v(0) * w(1)) - (v(1) * w(0)); 
        n(3) = n(0) * p1(0) + n(1) * p1(1) + n(2) * p1(2);

        Kp = n * n.transpose();
        Q += Kp;
      }
    }
    Qs.push_back(std::pair<Eigen::MatrixXd, int>(Q, i));
  }

  return Qs;
}

void updateViewer(igl::opengl::glfw::Viewer& viewer)
{
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
}

// Reset V and F to their original vertices and faces.
void reset(igl::opengl::glfw::Viewer& viewer)
{
  V = OV;
  F = OF;

  updateViewer(viewer);
}

int main(int argc, char *argv[])
{
  // Loads our Bunny mesh
  igl::readPLY("~/bunny10k.ply", OV, OF);
  igl::opengl::glfw::Viewer viewer;

  // Plot the mesh
  reset(viewer);

  // R/r - resets to the original untouched face
  // 1 - libigl's embedded decimation/edge collapse 
  // 2 - our implementation
  const auto &key_down =
  [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
  {
    switch(key)
    {
      case 'R':
      case 'r':
        reset(viewer);
        break;
      case '1':
        igl::decimate(V, F, F.size() * .20, V, F, j);
        updateViewer(viewer);
        break;
      case '2':
        igl::edges(F, E);
        computeOptimalContract(computeQ(V,F), E);
        // Colapse_Edge
      default: return false;
    }

    return true;
  };

  viewer.callback_key_down = key_down;
  viewer.launch();
}
