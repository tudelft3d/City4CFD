// Quadtree
//
// Copyright © 2022,
// Ken Arroyo Ohori    ken@ken.mx
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef Quadtree_h
#define Quadtree_h

template <class Kernel, class Point_cloud>
struct Quadtree_node {
  typename Kernel::FT x_min, x_max, y_min, y_max;
  unsigned int depth;
  std::vector<typename Point_cloud::Index> points;
  Quadtree_node *upper_left, *upper_right, *lower_left, *lower_right;
  
  Quadtree_node() {
    depth = 0;
    upper_left = NULL;
    upper_right = NULL;
    lower_left = NULL;
    lower_right = NULL;
  }
  
  Quadtree_node(unsigned int depth, typename Kernel::FT x_min, typename Kernel::FT x_max, typename Kernel::FT y_min, typename Kernel::FT y_max) {
    this->depth = depth;
    upper_left = NULL;
    upper_right = NULL;
    lower_left = NULL;
    lower_right = NULL;
    this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
  }
  
  void compute_extent(Point_cloud &point_cloud) {
    typename Kernel::FT x_min = point_cloud.point(*point_cloud.begin()).x();
    typename Kernel::FT x_max = point_cloud.point(*point_cloud.begin()).x();
    typename Kernel::FT y_min = point_cloud.point(*point_cloud.begin()).y();
    typename Kernel::FT y_max = point_cloud.point(*point_cloud.begin()).y();
    for (auto const &point: point_cloud.points()) {
      if (point.x() < x_min) x_min = point.x();
      if (point.x() > x_max) x_max = point.x();
      if (point.y() < y_min) y_min = point.y();
      if (point.y() > y_max) y_max = point.y();
    } this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
  }
  
  void insert_point(Point_cloud &point_cloud, typename Point_cloud::Index point) {
      points.push_back(point);
  }
  
  void optimise(Point_cloud &point_cloud, std::size_t bucket_size, unsigned int maximum_depth) {
    if (points.size() < bucket_size || depth+1 >= maximum_depth) return;
    typename Kernel::FT x_split = (x_min+x_max)/2.0;
    typename Kernel::FT y_split = (y_min+y_max)/2.0;
    upper_left = new Quadtree_node(depth+1, x_min, x_split, y_split, y_max);
    upper_right = new Quadtree_node(depth+1, x_split, x_max, y_split, y_max);
    lower_left = new Quadtree_node(depth+1, x_min, x_split, y_min, y_split);
    lower_right = new Quadtree_node(depth+1, x_split, x_max, y_min, y_split);
    
    for (auto point: points) {
      if (point_cloud.point(point).x() < x_split) {
        if (point_cloud.point(point).y() < y_split) lower_left->insert_point(point_cloud, point);
        else upper_left->insert_point(point_cloud, point);
      } else {
        if (point_cloud.point(point).y() < y_split) lower_right->insert_point(point_cloud, point);
        else upper_right->insert_point(point_cloud, point);
      }
    } points.clear();
    upper_left->optimise(point_cloud, bucket_size, maximum_depth);
    upper_right->optimise(point_cloud, bucket_size, maximum_depth);
    lower_left->optimise(point_cloud, bucket_size, maximum_depth);
    lower_right->optimise(point_cloud, bucket_size, maximum_depth);
  }
  
  bool intersects(typename Kernel::FT x_min, typename Kernel::FT x_max, typename Kernel::FT y_min, typename Kernel::FT y_max) {
//    std::cout << "\tBB X = [" << x_min << ", " << x_max << "] Y = [" << y_min << ", " << y_max << "] intersects X = [" << this->x_min << ", " << this->x_max << "] Y = [" << this->y_min << ", " << this->y_max << "]" << std::endl;
    if (x_max < this->x_min) return false;
    if (x_min > this->x_max) return false;
    if (y_min > this->y_max) return false;
    if (y_max < this->y_min) return false;
    return true;
  }
  
  void find_intersections(std::vector<Quadtree_node *> &nodes, typename Kernel::FT x_min, typename Kernel::FT x_max, typename Kernel::FT y_min, typename Kernel::FT y_max) {
    if (!points.empty()) nodes.push_back(this);
    if (upper_left != NULL && upper_left->intersects(x_min, x_max, y_min, y_max)) upper_left->find_intersections(nodes, x_min, x_max, y_min, y_max);
    if (upper_right != NULL && upper_right->intersects(x_min, x_max, y_min, y_max)) upper_right->find_intersections(nodes, x_min, x_max, y_min, y_max);
    if (lower_left != NULL && lower_left->intersects(x_min, x_max, y_min, y_max)) lower_left->find_intersections(nodes, x_min, x_max, y_min, y_max);
    if (lower_right != NULL && lower_right->intersects(x_min, x_max, y_min, y_max)) lower_right->find_intersections(nodes, x_min, x_max, y_min, y_max);
  }
  
  void print_info() const {
    std::cout << "Quadtree extent: X = [" << x_min << ", " << x_max << "] Y = [" << y_min << ", " << y_max << "]" << std::endl;
    std::vector<std::size_t> nodes_by_depth;
    Quadtree_node *biggest_node = this;
    std::list<Quadtree_node *> to_visit;
    to_visit.push_back(this);
    while (!to_visit.empty()) {
      if (to_visit.front()->points.size() > biggest_node->points.size()) biggest_node = to_visit.front();
      if (nodes_by_depth.size() < to_visit.front()->depth+1) nodes_by_depth.push_back(1);
      else nodes_by_depth[to_visit.front()->depth]++;
      if (to_visit.front()->upper_left != NULL) to_visit.push_back(to_visit.front()->upper_left);
      if (to_visit.front()->upper_right != NULL) to_visit.push_back(to_visit.front()->upper_right);
      if (to_visit.front()->lower_left != NULL) to_visit.push_back(to_visit.front()->upper_left);
      if (to_visit.front()->lower_right != NULL) to_visit.push_back(to_visit.front()->lower_right);
      to_visit.pop_front();
    } std::cout << "Quadtree nodes by depth:";
    for (auto const &n_nodes: nodes_by_depth) std::cout << " " << n_nodes;
    std::cout << std::endl;
    std::cout << "Biggest node: " << biggest_node->points.size() << " points at depth " << biggest_node->depth << " for X = [" << biggest_node->x_min << ", " << biggest_node->x_max << "] Y = [" << biggest_node->y_min << ", " << biggest_node->y_max << "]" << std::endl;
  }
  
  ~Quadtree_node() {
    if (upper_left != NULL) delete upper_left;
    if (upper_right != NULL) delete upper_right;
    if (lower_left != NULL) delete lower_left;
    if (lower_right != NULL) delete lower_right;
  }
};

#endif /* Quadtree_h */
