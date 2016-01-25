#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <vigra/impex.hxx>
#include <vigra/multi_array.hxx>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <just_watershed.hpp>

typedef boost::multi_array_types::index_range range;

int
num_files(boost::filesystem::path p) {

	int n = 0;
	for (boost::filesystem::directory_iterator i(p); i != boost::filesystem::directory_iterator(); i++)
		if (boost::filesystem::is_regular_file(*i))
			n++;

	return n;
}

std::vector<boost::filesystem::path>
files(boost::filesystem::path p) {

	std::vector<boost::filesystem::path> files;
	for (boost::filesystem::directory_iterator i(p); i != boost::filesystem::directory_iterator(); i++)
		if (boost::filesystem::is_regular_file(*i))
			files.push_back(*i);

	return files;
}

template <typename T>
void
read_slice_channel(
		T view,
		boost::filesystem::path aff_x_file) {

	vigra::ImageImportInfo info(aff_x_file.native().c_str());
	vigra::MultiArray<2, float> data(info.shape());
	vigra::importImage(aff_x_file.native().c_str(), data);

	for (int x = 0; x < info.shape()[0]; x++)
	for (int y = 0; y < info.shape()[1]; y++)
		view[x][y] = data(x, y);
}

template <typename T>
void
read_slice(
		T view,
		boost::filesystem::path aff_x_file,
		boost::filesystem::path aff_y_file,
		boost::filesystem::path aff_z_file) {

	read_slice_channel(
			view[ boost::indices[range()][range()][0] ],
			aff_x_file);
	read_slice_channel(
			view[ boost::indices[range()][range()][1] ],
			aff_y_file);
	read_slice_channel(
			view[ boost::indices[range()][range()][2] ],
			aff_z_file);
}

template <typename T>
void
write_slice(
		T view,
		std::string filename) {

	int size_x = view.shape()[0];
	int size_y = view.shape()[1];

	vigra::MultiArray<2, float> data(vigra::MultiArray<2, float>::difference_type(size_x, size_y));

	for (int x = 0; x < size_x; x++)
	for (int y = 0; y < size_y; y++)
		 data(x, y) = view[x][y];

	vigra::exportImage(
			data,
			filename.c_str());
}

struct dynamic_size_threshold {

	dynamic_size_threshold(int min_size, float threshold) :
		ms(min_size), t(threshold) {}

	std::size_t operator()(float v) const {

		// if not at least t affin, do not merge
		if (v < t) return 0;

		// otherwise, merge if smaller than ms
		return ms;
	}

	int ms;
	float t;
};

int main(int argc, char** argv) {

	if (argc < 8) {

		std::cout << "usage: watershed <aff_x_dir> <aff_y_dir> <aff_z_dir> <t_l> <t_h> <t_s> <ms>" << std::endl;
		return 1;
	}

	std::string aff_x_dir = argv[1];
	std::string aff_y_dir = argv[2];
	std::string aff_z_dir = argv[3];
	float t_l = boost::lexical_cast<float>(argv[4]);
	float t_h = boost::lexical_cast<float>(argv[5]);
	float t_s = boost::lexical_cast<float>(argv[6]);
	int   ms  = boost::lexical_cast<float>(argv[7]);

	std::cout
			<< "Performing affinity graph watershed on volumes "
			<< aff_x_dir << ", " << aff_y_dir << ", " << aff_z_dir
			<< std::endl;

	boost::filesystem::path aff_x_path(aff_x_dir);
	boost::filesystem::path aff_y_path(aff_y_dir);
	boost::filesystem::path aff_z_path(aff_z_dir);

	int size_z = num_files(aff_x_path);
	if (size_z != num_files(aff_y_dir) || size_z != num_files(aff_z_dir)) {

		std::cerr << "directories contain different number of files" << std::endl;
		return 1;
	}

	if (size_z == 0) {

		std::cerr << "directories contain no files" << std::endl;
		return 1;
	}

	std::vector<boost::filesystem::path> aff_x_files = files(aff_x_path);
	std::vector<boost::filesystem::path> aff_y_files = files(aff_y_path);
	std::vector<boost::filesystem::path> aff_z_files = files(aff_z_path);
	aff_x_files.resize(1); // one section only
	aff_y_files.resize(1); // one section only
	aff_z_files.resize(1); // one section only
	size_z = 1;

	std::sort(aff_x_files.begin(), aff_x_files.end());
	std::sort(aff_y_files.begin(), aff_y_files.end());
	std::sort(aff_z_files.begin(), aff_z_files.end());

	vigra::ImageImportInfo info(aff_x_files[0].native().c_str());
	int size_x = info.width();
	int size_y = info.height();

	std::cout << "reading affinity graph of size " << size_x << "x" << size_y << "x" << size_z << std::endl;

	affinity_graph_ptr<float> aff(
			new affinity_graph<float>(
					boost::extents[size_x][size_y][size_z][3],
					boost::fortran_storage_order()));

	for (int z = 0; z < size_z; z++) {

		auto slice = (*aff)[ boost::indices[range()][range()][z][range()] ];
		read_slice(slice, aff_x_files[z], aff_y_files[z], aff_z_files[z]);
	}

	std::cout << "performing simple_watershed" << std::endl;

	std::vector<std::size_t> counts;
	auto result = simple_watershed<uint32_t>(aff, t_l, t_h, counts);

	auto segmentation = result.first;
	int num_segments  = result.second;

	std::cout << "found " << num_segments << " segments" << std::endl;

	auto rg = get_region_graph<uint32_t,float>(aff, segmentation, num_segments);

	std::cout << "performing region merging" << std::endl;

	// I guess the last parameter is to discard regions smaller than that
	//merge_segments_with_function(result.first, rg, counts, 
	//dynamic_size_threshold(t_s, ms), 0);

	for (int z = 0; z < size_z; z++) {

		auto slice = (*segmentation)[ boost::indices[range()][range()][z] ];
		std::stringstream filename;
		filename << "watershed_" << std::setw(5) << std::setfill('0') << z << "_" << t_l << "_" << t_h << "_" << t_s << "_" << ms << ".tif";
		write_slice(slice, filename.str());
	}
}
