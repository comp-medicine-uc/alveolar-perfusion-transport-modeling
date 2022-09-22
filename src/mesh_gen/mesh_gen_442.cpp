#include "../../../cgal/CGAL/Kernel_23/include/Exact_predicates_inexact_constructions_kernel.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/Mesh_triangulation_3.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/Mesh_complex_3_in_triangulation_3.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/Mesh_criteria_3.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/Labeled_mesh_domain_3.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/make_mesh_3.h"
#include "../../../cgal/CGAL/Kernel_23/include/CGAL/Image_3.h"
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
int main(int argc, char* argv[])
{
    const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("rve.inr");
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }
    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);
    // Mesh criteria
    Mesh_criteria criteria(facet_angle = 30, facet_size = 6, facet_distance = 4,
        cell_radius_edge_ratio = 3, cell_size = 8);
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    // Output
    std::ofstream medit_file("out.mesh");
    c3t3.output_to_medit(medit_file);
    return 0;
}