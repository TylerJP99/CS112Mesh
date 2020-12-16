
=======
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
<<<<<<< HEAD
>>>>>>> 899536b38edb22e24e0cd0176bc56738111f37ea
=======
#include <igl/decimate.h>
>>>>>>> 7b0ad1c591be9f4e7b5f934120f3138f023bdbc8

Eigen::MatrixXd V, OV;
Eigen::MatrixXi F, OF;
Eigen::VectorXi j;

// Reset V and F to their original vertices and faces.
auto reset(igl::opengl::glfw::Viewer& viewer)
{
  V = OV;
  F = OF;

  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
}

// Maybe have a button that'll simplify the mesh by a certain amount?
// And a reset button to restore the mesh to the original unsimplified. My idea would be :
// r - resets to the original untouched face
// 1 - libigl's embedded decimation/edge collapse
// 2 - our implementation
int main(int argc, char *argv[])
{
  // Loads our Bunny mesh
  igl::readPLY("C:/Users/Donny/OneDrive/Documents/GitHub/personal-example/bunny10k.ply", OV, OF);
  igl::opengl::glfw::Viewer viewer;

  // Plot the mesh
  reset(viewer);

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
        // Make the number of faces percentage based not hardcoded at 1000
        igl::decimate(OV, OF, 1000, V, F, j);
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.data().set_face_based(true);
        break;
      default: return false;
    }

    return true;
  };


  viewer.callback_key_down = key_down;
  // viewer.data().set_mesh(V, F);
  // viewer.data().set_face_based(true);
  viewer.launch();
}
