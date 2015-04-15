/* Definitions of major classes */ 

#include "rtnode.h"

class RTree {
	public:
		RTree(int entry_num);//by default, dimension is 2
		RTree(int entry_num, int dim);
		~RTree();

	private:
		bool same_entry(const Entry& e1, const Entry& e2);
		bool overlap(const BoundingBox box1, const BoundingBox box2);
		void update_mbr(BoundingBox& mbr, const BoundingBox& new_mbr);
		BoundingBox get_mbr(Entry* entry_list, int len);
		int area(const BoundingBox& mbr);
		void swap_entry(Entry* entry_list, int id1, int id2);
		int area_inc(const BoundingBox& mbr, const BoundingBox& entry_mbr);
		void linear_pick_seeds(Entry* entry_list, int len, int& m1, int& m2);
		RTNode* find_leaf(RTNode* node, RTNode** stack, int* entry_idx, int& stack_size, const Entry& record);
		RTNode* choose_leaf(RTNode** stack, int* entry_idx, int& stack_size, const Entry& record, int dest_level);
		void adjust_tree(RTNode** stack, int* entry_idx, int size);
		void query_range(const RTNode* node, const BoundingBox mbr, int& result_cnt, int& node_travelled);
		bool query_point(const RTNode* node, const BoundingBox& mbr, Entry& result);
		bool insert(const Entry& e, int dest_level);
		void stat(RTNode* node, int& record_cnt, int& node_cnt);
		void print_node(RTNode* node, int indent_level);

	public:
		void stat();
		void print_tree();
		bool insert(const vector<int>& coordinate, int rid);
		void query_range(const BoundingBox& mbr, int& result_count, int& node_travelled);
		bool query_point(const vector<int>& coordinate, Entry& result);
		bool tie_breaking(const BoundingBox& box1, const BoundingBox& box2);
		bool del(const vector<int>& coordinate);

	private:
		int max_entry_num;
		int dimension;
		RTNode* root;
};
