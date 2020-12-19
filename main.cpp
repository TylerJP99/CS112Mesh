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
std::vector<std::pair<Eigen::Vector4d, int>> Vhats;
Eigen::Matrix<int, Eigen::Dynamic, 2> E;

void edge_Collapse(
  Eigen::MatrixXd &V, Eigen::MatrixXi &F,
  Eigen::Vector4d Vh, std::pair<int, int> e)
  {

    Eigen::MatrixXd NV;
    Eigen::MatrixXi NF;
    int r = 0;
    int vi;
    bool newV_added = false;
    Eigen::Vector3d vh = Vh.head<3>();
    NV.resize(V.rows()-1, 3);
    NF.resize(F.rows()-2, 3);

    for (unsigned int i = 0; i < V.rows(); i++){
      if (i==e.first || i==e.second){
        if (newV_added){

        }
        else{
          newV_added = true;
          NV.row(r) = vh;
          vi = r;
          r++;
        }
      }
      else{
        NV.row(r) = V.row(i);
        r++;
      }
    }
    r = 0;

    for (unsigned int i = 0; i < F.rows(); i++){
      if ((F(i, 0)==e.first && F(i, 1)==e.second)
      || (F(i, 0)==e.second && F(i, 1)==e.first)
      || (F(i, 0)==e.first && F(i, 2)==e.second)
      || (F(i, 0)==e.second && F(i, 2)==e.first)
      || (F(i, 1)==e.first && F(i, 2)==e.second)
      || (F(i, 1)==e.second && F(i, 2)==e.first))
      {

      }
      else if ((F(i, 0)==e.first || F(i, 0)==e.second)
      || (F(i, 1)==e.first || F(i, 1)==e.second)
      || (F(i, 2)==e.first || F(i, 2)==e.second)){
        if ((F(i, 0)==e.first || F(i, 0)==e.second)){
          NF.row(r) = F.row(i);
          NF(r, 0) = vi;
          r++;
        }
        else if ((F(i, 1)==e.first || F(i, 1)==e.second)){
          NF.row(r) = F.row(i);
          NF(r, 1) = vi;
          r++;
        }
        else if ((F(i, 2)==e.first || F(i, 2)==e.second)){
          NF.row(r) = F.row(i);
          NF(r, 2) = vi;
          r++;
        }
      }
      else{
        NF.row(r) = F.row(i);
        r++;
      }
    }
    V = NV;
    F = NF;
  }


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
    double cost = Vhat.transpose() * Qhat * Vhat;

    Q.emplace(cost, i, 0);
    std::pair<Eigen::Vector4d, int> vh;
    vh.first = Vhat;
    vh.second = i;
    Vhats.push_back(vh);
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

void Decimate(Eigen::MatrixXd &V, Eigen::MatrixXi &F, int m_f){
  std::tuple<double, int , int> p;
  std::pair<Eigen::Vector4d, int> Vh;
  std::pair<int, int> e;
  Eigen::Vector4d vh;
  int m = 0;
  while(true){
    if (F.rows() <= m_f){
      break;
    }
    else if (m==10000){
      break;
    }
    else{
      igl::edges(F, E);
      computeOptimalContract(computeQ(V,F), E);
      if (Q.size()==0){
        break;
      }
      p = Q.top();
      Q.pop();
      for (unsigned int i = 0; i<Vhats.size(); i++){
        if (Vhats[i].second == std::get<1>(p)){
          Vh = Vhats[i];
          break;
        }
      }
      e.first = E(std::get<1>(p), 0);
      e.second = E(std::get<1>(p), 1);
      vh = Vh.first;
      edge_Collapse(V,F,vh,e);
    }
  m++;
  }
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
        igl::decimate(V, F, F.rows() * .20, V, F, j);
        updateViewer(viewer);
        break;
      case '2':
        Decimate(V, F, F.rows() * .20);
        updateViewer(viewer);
        // Colapse_Edge
      default: return false;
    }

    return true;
  };

  viewer.callback_key_down = key_down;
  viewer.launch();
}
