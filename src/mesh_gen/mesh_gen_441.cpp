#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <functional>
typedef float Image_word_type;
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
int main(int argc, char* argv[])
{
    const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/skull_2.9.inr");
    // Load image
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }
    Mesh_domain domain =
        Mesh_domain::create_gray_image_mesh_domain(image, 2.9f, 0.f);
    // Mesh criteria
    Mesh_criteria criteria(facet_angle = 30, facet_size = 6, facet_distance = 2,
        cell_radius_edge_ratio = 3, cell_size = 8);
    // Meshing
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    // Output
    std::ofstream medit_file("out.mesh");
    c3t3.output_to_medit(medit_file);
    return 0;
}