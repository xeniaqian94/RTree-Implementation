
#include "rtnode.h"

class RTree {
	public:
		RTree(int entry_num);//by default, dimension is 2
		RTree(int entry_num, int dim);
		~RTree();

	private:
		//DO NOT change below three functions.
		//compute the MBR of all the entries in entry list
		BoundingBox get_mbr(Entry* entry_list, int len);
		//show statistics, called by public stat()
		void stat(RTNode* node, int& record_cnt, int& node_cnt);
		//print a node
		void print_node(RTNode* node, int indent_level);

	public:
		//DO NOT change below three functions
		void stat();//print statistics of the R-tree
		void print_tree();//print the R-tree
		bool tie_breaking(const BoundingBox& box1, const BoundingBox& box2);
		//PLEASE implement below three functions
		bool insert(const vector<int>& coordinate, int rid);
    
        RTNode* choose_leaf(const BoundingBox& point_box, int rid);
        void split_node(RTNode*& leaf_l, RTNode*& leaf_ll,Entry* entry_list);
        void linear_pick_seed(Entry* entry_list,int& seed_x, int& seed_y);
        
        void find_pairs(int dimension, Entry* entry_list, double& max, int& x, int &y);
        void adjust_tree(RTNode*& leaf_l, RTNode*& leaf_ll,RTNode*& p, RTNode*& pp);
    
		void query_range(const BoundingBox& mbr, int& result_count, int& node_travelled);
        void query_range_subtree(RTNode* root, const BoundingBox& mbr, int& result_count, int& node_travelled);
    
		bool query_point(const vector<int>& coordinate, Entry& result);
        bool query_point_subtree(RTNode* root, const vector<int>& coordinate, Entry& result);
    
	private:
		int max_entry_num;
		int dimension;
		RTNode* root;
};
