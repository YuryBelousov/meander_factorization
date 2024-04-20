#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>

/*
This file contains the class which is used to enumerate irreducible meanders 
(after minor corrections it could be used to directly enumerate all meanders).


							ALGORITHM DESCRIPTION
We use the terminology from the corresponding paper. We assume that there are n points on l through which 
m can pass through. The algorithm recursively visits all available points, checking that it can end up with an 
almost irreducible meander. To check the available points, we use a graph whose vertices are the regions 
into which l and m divide the disk, and whose edges are the unvisited points on l. At each step we only need 
to check that:
(1) if we remove multiple edges and vertices of degree 0, the resulting graph is connected;
(2) the sequence of visited points corresponds to an almost irreducible meander.

To speed up the calculation, we omit some checks that are generally required.
*/

constexpr int meanless_value = -1;	//Special value, used in several situations
constexpr int max_n = 24;			//Maximum number of intersections

class submeander {
private:
	//Meander code
	int code[max_n];	//Permutation of the meander
	int code_size;		//Number of added intersections

	int tail_index;		//Label of the last visited area
	int n;				//Maximum number in the code

	bool next_point_odd;			//Parity of the next point
	bool not_visited[max_n];		//Shows whether a particular point was visited

	int upper_areas[max_n];			//Label of the area above intersection point
	int lower_areas[max_n];			//Label of the area below intersection point
	int areas_hight[max_n + 3];		//Height of the area with corresponding label

	int edges[max_n / 2 + 2][max_n / 2 + 2];	//edges[i][j] is the number of edges from upper_areas[i]/2 to (lower_areas[j]-1)/2
	int verticies[max_n + 3];					//Degrees of the areas (interpret them as vertecies in the dual graph)
	int edge_num;								//Total number of different edges
	int vert_num;								//Total number of verticies with non zero incidence index

	mutable int avaliable_points[max_n / 2 + 1];	//The array stores labels of vertecies that we can visit in the next step

private:
	void remove_edge(int from, int to) noexcept {
		if (--edges[from / 2][(to - 1) / 2] == 0) {
			--edge_num;
			if (--verticies[from] == 0)
				--vert_num;
			if (--verticies[to] == 0)
				--vert_num;
		}
	}

	void add_edge(int from, int to) noexcept {
		if (++edges[from / 2][(to - 1) / 2] == 1) {
			++edge_num;
			if (++verticies[from] == 1)
				++vert_num;
			if (++verticies[to] == 1)
				++vert_num;
		}
	}

public:
	//Create an empty submeander.
	submeander() noexcept
		: n(0), next_point_odd(false), tail_index(meanless_value), code_size(0), edge_num(1), vert_num(2) {
		//Create upper and lower half-plane
		for (int i = 0; i < max_n; ++i) {
			not_visited[i] = true;
			lower_areas[i] = 1;
		}
		std::memset(&edges, 0, sizeof(edges));
		std::memset(&verticies, 0, sizeof(verticies));
		std::memset(&upper_areas, 0, sizeof(upper_areas));

		areas_hight[0] = areas_hight[1] = max_n + 3;
		verticies[upper_areas[0]] = verticies[lower_areas[0]] = 1;
		edges[0][0] = max_n + 2;
	}

	//Create submeander based on a given permutation. Note: the labels of the intersection starts with zero.
	submeander(const std::vector<int>& perm) noexcept : submeander() {
		for (const auto& pt : perm)
			go_throught_point(pt);
	}

	//Modify submeander by passing through new point.
	//Note: there is no check that this operation is correct.
	void go_throught_point(int new_point) noexcept {
		n = std::max(n, new_point + 1);
		//Choose half plane that will change
		int* half_plane = next_point_odd ? lower_areas : upper_areas;
		//Add new area
		int new_area_height = std::abs(tail_index - new_point);
		areas_hight[code_size + 2] = new_area_height;
		//Remove an edge in new_point
		remove_edge(upper_areas[new_point], lower_areas[new_point]);
		//Points between now lead to new area
		for (int i = std::min(tail_index, new_point) + 1; i < std::max(tail_index, new_point); ++i)
			if (areas_hight[half_plane[i]] > new_area_height) {
				if (not_visited[i]) {
					remove_edge(upper_areas[i], lower_areas[i]);
					if (next_point_odd)
						add_edge(upper_areas[i], code_size + 2);
					else
						add_edge(code_size + 2, lower_areas[i]);
				}
				half_plane[i] = code_size + 2;
			}
		//Delete new_point from not_visited 
		not_visited[new_point] = false;
		//Change other info
		tail_index = new_point;
		next_point_odd = !next_point_odd;
		code[code_size] = new_point;
		++code_size;
	}

	//Delete the last visited point from the submeander. 
	//Note: to use it to enumerate all meanders, you should uncomment additional check.
	void step_back() noexcept {
		// Note: to speed up calculations, this function only works when code_size > 2.
		/*if (code_size < 2)
			tail_index = -1;
		else*/
			tail_index = code[code_size - 2];

		if (n == code[code_size - 1] + 1) {
			n = meanless_value;
			for (int i = 0; i < code_size - 1; ++i)
				n = std::max(n, code[i] + 1);
		}
		not_visited[code[code_size - 1]] = true;
		int* half_plane = next_point_odd ? upper_areas : lower_areas;
		int prev_area = half_plane[code[code_size - 1]];
		for (int i = std::min(tail_index, code[code_size - 1]); i < std::max(tail_index, code[code_size - 1]); ++i)
			if (half_plane[i] == code_size + 1) {
				half_plane[i] = prev_area;
				if (not_visited[i]) {
					add_edge(upper_areas[i], lower_areas[i]);
					if (next_point_odd)
						remove_edge(code_size + 1, lower_areas[i]);
					else
						remove_edge(upper_areas[i], code_size + 1);
				}
			}
		next_point_odd = !next_point_odd;
		--code_size;
		add_edge(upper_areas[code[code_size]], lower_areas[code[code_size]]);
	}

	//Find all points available for the go_through function.
	//Note: to use it to enumerate all meanders, you should uncomment additional check.
	const int* const get_avaliable_points(/*bool irreducibility = true*/) const noexcept {
		if (vert_num != edge_num + 1) {
			avaliable_points[0] = meanless_value;
			return avaliable_points;
		}
		//Choose half plane where we a looking for avaliable points
		const int* const half_plane = next_point_odd ? lower_areas : upper_areas;
		//Number of tail area
		int tail_area = (tail_index == meanless_value) ? 0 : half_plane[tail_index];
		int ind = -1;
		for (int i = next_point_odd ? 1 : 0; i < max_n; i += 2)
			if (not_visited[i] &&
				half_plane[i] == tail_area &&
				(/*!irreducibility ||*/ irreducibility_save(i)))
				avaliable_points[++ind] = i;
		avaliable_points[++ind] = meanless_value;
		return avaliable_points;
	}

	//This function checks that adding a new point doesn't break the irreducibility.
	//Note: to use it to enumerate all meanders, you should uncomment additional check.
	bool irreducibility_save(int pt) const noexcept {
		if (code_size < 2)
			return true;

		//This check is useful in the general case, 
		//but we only consider meanders that do not start at the first point, so we omit it to speed up the calculations.
		/*if (code[0] == 0)
			return false;*/
		
		//We start at the end of the code and try to find sequence{ a_i } of length k(k > 2) 
		// with the property that max {a_i} - min{ a_i } = k.
		int new_n = std::max(n, pt + 1),
			min_p = pt,
			max_p = pt;
		for (int i = code_size - 1; i > -1; --i) {
			if (min_p > code[i])
				min_p = code[i];
			else if (max_p < code[i])
				max_p = code[i];
			if (code_size - i < 2)
				continue;
			if ((i == 0) && (code_size + 1 == new_n))
				break;
			if ((code_size - i) == max_p - min_p)
				return false;
		}
		return true;
	}

	//Check that it is a meander (i.e. all the first points have been visited and 
	// the curve can be continued to the boundary of the disc without adding new intersections).
	bool finished() const noexcept {
		if (code_size < n)
			return false;
		if (next_point_odd)
			return (lower_areas[tail_index] == 1);
		else
			return (upper_areas[tail_index] == 0);
	}

	int get_n() const noexcept {
		return n;
	}

	int get_cups() const noexcept {
		int cups_number = 0;
		for (int i = 1; i < code_size; ++i)
			if (std::abs(code[i] - code[i - 1]) == 1)
				++cups_number;
		return cups_number;
	}

	int get_first_point() const noexcept {
		return code[0];
	}

	int get_last_point() const noexcept {
		return code[code_size - 1];
	}

	int get_length() const noexcept {
		int length = 0;
		for (int i = 1; i < code_size; ++i)
			length += std::abs(code[i] - code[i - 1]);
		return length;
	}

	//Comparison operator required to store submeanders in set.
	friend bool operator<(const submeander& m1, const submeander& m2) {
		if (m1.code_size != m2.code_size)
			return m1.code_size < m2.code_size;
		for (size_t i = 0; i < m1.code_size; ++i)
			if (m1.code[i] != m2.code[i])
				return m1.code[i] < m2.code[i];
		return false;
	}

	//Print submeander permutation in a stream.
	friend std::ostream& operator<<(std::ostream& out, const submeander& m) noexcept {
		for (size_t i = 0; i < std::min(m.code_size, max_n); ++i)
			out << m.code[i] + 1 << " ";
		return out;
	}

	void make_tikz_picture(std::ostream& out, double scale = 3, bool draw_circle = true) const noexcept {
		double right_border = 1;
		int max_h_up = 0,
			max_h_bottom = 0;
		//Find y-coordinate of highest and lowest arcs
		for (int i = 0; i < code_size - 1; ++i) {
			if (i % 2 == 1)
				max_h_up = std::max(max_h_up, abs(code[i + 1] - code[i]));
			else
				max_h_bottom = std::max(max_h_bottom, abs(code[i + 1] - code[i]));
		}
		double
			curvature = 12.56637,
			h = right_border / (double)(code_size + 1),
			circle_center = right_border / 2.,
			line_x_start = 0,
			line_x_end = right_border,
			circle_radius = (line_x_end - line_x_start) / 2.,
			m_y_start = h * (max_h_up + 2.) / 3,
			m_y_end = (code_size % 2 == 0) ? h * (max_h_up + 2) / 3. : -h * (max_h_bottom + 2) / 3.,
			m_x_start = (draw_circle) ? -std::sqrt(circle_radius * circle_radius - m_y_start * m_y_start) + circle_center : line_x_start,
			m_x_end = (draw_circle) ? std::sqrt(circle_radius * circle_radius - m_y_end * m_y_end) + circle_center : line_x_end;
		out << "\\begin{tikzpicture}[scale = " << scale << "]\n";
		//Draw horizontal line
		out << "\\draw[thick] (" << line_x_start << ", 0) to (" << line_x_end << ", 0);\n";
		//Draw curve 
		out << "\\draw[ultra thick] (" << m_x_start << ", " << m_y_start << ")";
		int angle = 90;
		for (int i = 0; i < code_size; ++i) {
			out << " to[out = " << ((i == 0) ? 0 : angle) << ", in = " << angle
				<< ", distance = " << h * abs(code[i] - ((i == 0) ? -1 : code[i - 1])) * curvature << "] ("
				<< h * (code[i] + 1) << ", 0)\n";
			angle *= -1;
		}
		out << "to[out = " << angle << ", in = 180, distance = " <<
			h * abs(code_size - code[code_size - 1]) * curvature << "] (" <<
			m_x_end << ", " << m_y_end << ");\n";
		if (draw_circle) {
			double point_radius = sqrt(3. / scale) / 60.;
			out << "\\draw[help lines] (" << circle_center << ", 0) circle (" << circle_radius << ");\n";
			out << "\\draw[fill] (" << m_x_start << ", " << m_y_start << ") circle (" << (point_radius) << ");\n";
			out << "\\draw[fill] (" << m_x_end << ", " << m_y_end << ") circle (" << (point_radius) << ");\n";
			out << "\\draw[fill] (" << line_x_start << ", " << 0 << ") circle (" << (point_radius) << ");\n";
			out << "\\draw[fill] (" << line_x_end << ", " << 0 << ") circle (" << (point_radius) << ");\n";
		}
		out << "\\end{tikzpicture}\n";
	}
};