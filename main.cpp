
=======
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
>>>>>>> 899536b38edb22e24e0cd0176bc56738111f37ea

int main(int argc, char *argv[])
{
  // Loads our Bunny mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readPLY("C:/Users/Donny/OneDrive/Documents/GitHub/personal-example/bunny10k.ply", V, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
