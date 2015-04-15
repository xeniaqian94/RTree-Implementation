/* Implementations of R tree */
#include <cmath>
#include "rtree.h"


const double EPSILON = 1E-10;

RTree::RTree(int entry_num)
{
	max_entry_num = entry_num;
	dimension = 2;//by default
	root = new RTNode(0, entry_num);
}

RTree::RTree(int entry_num, int dim)
{
	max_entry_num = entry_num;
	dimension = dim;//by default
	root = new RTNode(0, entry_num);
}

RTree::~RTree()
{
	delete root;
	root = NULL;
}


//
// Check whether two entries are the same.
// Return true if same, otherwise false.
//
bool RTree::same_entry(const Entry& e1, const Entry& e2)
{
	const BoundingBox& mbr1 = e1.get_mbr();
	const BoundingBox& mbr2 = e2.get_mbr();
	
	return mbr1.is_equal(mbr2);
}


//
// Check whether two boundingboxs overlap.
// Return true if so, otherwise false.
//
bool RTree::overlap(const BoundingBox box1, const BoundingBox box2)
{
	return box1.is_intersected(box2);
}


//
// Update the current MBR ``mbr'' with the merging result of ``mbr'' and ``new_mbr''
//
void RTree::update_mbr(BoundingBox& mbr, const BoundingBox& new_mbr)
{
	mbr.group_with(new_mbr);
}


//
// Calcuate the MBR of a set of entries, of size ``len''.

BoundingBox RTree::get_mbr(Entry* entry_list, int len)
{
	BoundingBox mbr(entry_list[0].get_mbr());
	for (int i = 1; i < len; i++) {        
		mbr.group_with(entry_list[i].get_mbr());
	}
	return mbr;
}


//
// Return the area of a boundingbox ``mbr''.
//
int RTree::area(const BoundingBox& mbr)
{
	return mbr.get_area();
}


//
// Swap two entries: entry_list[id1] and entry_list[id2].
//
void RTree::swap_entry(Entry* entry_list, int id1, int id2)
{
	Entry temp = entry_list[id1];
	entry_list[id1] = entry_list[id2];
	entry_list[id2] = temp;
}


//
// Calculate the area enlarged by add the new entry to the existing MBR.
//
int RTree::area_inc(const BoundingBox& mbr, const BoundingBox& entry_mbr)
{
	BoundingBox new_mbr(mbr);
	new_mbr.group_with(entry_mbr);
	return area(new_mbr) - area(mbr);
}

//
// Linear Pick Seeds algorithm for Lienar Cost Algorithm.
//
void RTree::linear_pick_seeds(Entry* entry_list, int len, int& m1, int& m2)
{
	int iminx = 0, imaxx = 0, iminy = 0, imaxy = 0;    
	int dim = entry_list[0].get_mbr().get_dim();

	//extreme pairs for each dimension
	//for a pair, first element is the entry with highest low side, second is the entry with lowest high side
	vector<pair<int,int> > extremePairs;
	//initialize entreme pairs
	for (int i = 0; i < dim; i++)
	{
		pair<int, int> extremePair(0, 0);//the first entry by default
		extremePairs.push_back(extremePair);
	}
	// pick the two entries with the largest gap on each dimension
	//for every entry
	for (int i = 1; i < len; i++) {
		//for every dimension
		for (int j = 0; j < dim; j++)
		{
			pair<int, int> extremePair = extremePairs[j];
			BoundingBox ithMBR = entry_list[i].get_mbr();

			//get highest low side on j-th dimension
			//the MBR of entry that has highest low side on dimension j
			BoundingBox highestLowEntryMBR = entry_list[extremePair.first].get_mbr();
			if (ithMBR.get_lowestValue_at(j) > highestLowEntryMBR.get_lowestValue_at(j)) {
				extremePairs[j].first = i;
			}
			else if (ithMBR.get_lowestValue_at(j) == highestLowEntryMBR.get_lowestValue_at(j)) {
				if (tie_breaking(ithMBR, highestLowEntryMBR)) {
					extremePairs[j].first = i;
				}
			}

			//get lowest high side on j-th dimension
			BoundingBox lowestHighEntryMBR = entry_list[extremePair.second].get_mbr();
			if (ithMBR.get_highestValue_at(j) < lowestHighEntryMBR.get_highestValue_at(j))
			{
				extremePairs[j].second = i;
			}
			else if (ithMBR.get_highestValue_at(j) == lowestHighEntryMBR.get_highestValue_at(j))
			{
				if (tie_breaking(ithMBR, lowestHighEntryMBR)) {
					extremePairs[j].second = i;
				}
			}
		}

	}
	BoundingBox box = get_mbr(entry_list, len);
	
	//for each dimension, find the greatest normalized separation and store the respective pair in m1 and m2
	//init
	double greatestNormalizedSeparation = -1; // the normal value of this should be >= 0
	m1 = -1;
	m2 = -1;
	for (int j = 0; j < dim; j++)
	{
		double normalizedJdimSeparation = 0;
		double delta = box.get_highestValue_at(j) - box.get_lowestValue_at(j);

		pair<int, int> extremePair = extremePairs[j];
		if (delta != 0)
		{
			normalizedJdimSeparation = 
				abs(entry_list[extremePair.first].get_mbr().get_lowestValue_at(j) 
					- entry_list[extremePair.second].get_mbr().get_highestValue_at(j)) * 1.0 
				/ delta;
		}
		if (greatestNormalizedSeparation - normalizedJdimSeparation >= -EPSILON)
		{
		}
		else {
			m1 = extremePair.first;
			m2 = extremePair.second;
			greatestNormalizedSeparation = normalizedJdimSeparation;
		}
	}
	//tie breaking
	if(m1 == m2) {
		m2 = (m1 == 0 ? 1 : 0);
		for (int i = 1; i < len; i++) {
			if(i != m1 && i != m2) {
				if(tie_breaking(entry_list[i].get_mbr(), entry_list[m2].get_mbr()))
					m2 = i;
			}
		}
	}
}


//
// Find the leaf node and delete the ``record''.
//
RTNode* RTree::find_leaf(RTNode* node, RTNode** stack, int* entry_idx, int& stack_size, const Entry& record)
{
	if (node->level == 0) {
		for (int i = 0; i < node->entry_num; i++) {
			if (overlap(node->entries[i].get_mbr(), record.get_mbr())) {
				swap_entry(node->entries, i, node->entry_num-1); // move the record the the end to indicate ``deleted''
				node->entry_num--;
				return node;
			}
		}
	}
	else {
		for (int i = 0; i < node->entry_num; i++) {
			if (overlap(node->entries[i].get_mbr(), record.get_mbr())) {
				stack[stack_size] = node;
				entry_idx[stack_size] = i;
				stack_size++;
				RTNode* ret = find_leaf(node->entries[i].get_ptr(), stack, entry_idx, stack_size, record);
				if (ret != NULL) {
					return ret;
				}
				stack_size--;
			}
		}
	}
	return NULL;
}

//
// Find the node to insert the new entry ``e'' at the specified level ``dest_level''.
// In particular, find the leaf node for new record if ``dest_level == 0''.
//
RTNode* RTree::choose_leaf(RTNode** stack, int* entry_idx, int& stack_size, const Entry& e, int dest_level)
{
	RTNode* node = root;
	while (node->level != dest_level) {
		int min_idx = 0;
		int min_enlargement = area_inc(node->entries[0].get_mbr(), e.get_mbr());
		for (int i = 1; i < node->entry_num; i++) {
			// compare with other entries
			int cur_enlargement = area_inc(node->entries[i].get_mbr(), e.get_mbr());
			if (cur_enlargement < min_enlargement) {
				min_idx = i;
				min_enlargement = cur_enlargement;
			}
			else if (cur_enlargement == min_enlargement) {
				// do not need to change min_enlargement as they are the same.
				int cur_area = area(node->entries[i].get_mbr());
				int min_area = area(node->entries[min_idx].get_mbr());
				// select the one with min area.
				if (cur_area < min_area) {
					min_idx = i;
				}
				else if (cur_area == min_area) {
					// tie breaking
					if (tie_breaking(node->entries[i].get_mbr(), node->entries[min_idx].get_mbr())) {
						min_idx = i;
					}
				}
			}
		}
		//this->print_node(node, 4);
		stack[stack_size] = node;
		entry_idx[stack_size] = min_idx;
		stack_size++;
		node = node->entries[min_idx].get_ptr();
	}
	return node;
}


//
// Adjust the MBR of nodes involved in insertion.
//
void RTree::adjust_tree(RTNode** stack, int* entry_idx, int size)
{
	while (size > 0) {
		size--;
		RTNode* node = stack[size]->entries[entry_idx[size]].get_ptr();
		
		stack[size]->entries[entry_idx[size]].set_mbr(get_mbr(node->entries, node->entry_num));
	}
}


//
// Helper function for query_range(), with range specified in ``mbr''.
// Return: number of results in ``result_cnt''.
//		number of R-tree nodes traveled in ``node_traveled''.
void RTree::query_range(const RTNode* node, const BoundingBox mbr, int& result_cnt, int& node_traveled)
{
	node_traveled++;
	if (node->level == 0) {
		for (int i = 0;i < node->entry_num;i++) {
			if (overlap(node->entries[i].get_mbr(), mbr)) {
				result_cnt++;
			}
		}
	} else {
		for (int i = 0;i < node->entry_num; i++) {
			if (overlap(node->entries[i].get_mbr(), mbr)) {
				query_range(node->entries[i].get_ptr(), mbr, result_cnt, node_traveled);
			}
		}
	}
}


//
// Helper function for point_query().
//
bool RTree::query_point(const RTNode* node, const BoundingBox& mbr, Entry& result)
{
	if (node->level == 0) {
		for (int i = 0; i < node->entry_num; i++) {
			if (overlap(node->entries[i].get_mbr(), mbr)) {
				result = node->entries[i];
				return true;
			}
		}
	}
	else {
		for (int i = 0; i < node->entry_num; i++) {
			if (overlap(node->entries[i].get_mbr(), mbr)) {
				if (query_point(node->entries[i].get_ptr(), mbr, result)) {
					return true;
				}
			}
		}
	}
	return false;
}	


bool RTree::insert(const vector<int>& coordinate, int rid)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}
	//a point is also modeled by a mbr.
	BoundingBox mbr(coordinate, coordinate);
	Entry e(mbr, rid);
	return insert(e, 0);
}


//
// Helper function for insertion.
//
bool RTree::insert(const Entry& e, int dest_level)
{
	Entry dummy;
	if (dest_level == 0 && query_point(root, e.get_mbr(), dummy)) {
		return false; 
	} 

	// stack contains the path to the leaf (not including the leaf node).
	RTNode** stack = new RTNode*[root->level];
	int stack_size = 0;
	// entry_idx contains the index of each entry in the node from the path.
	int* entry_idx = new int[root->level];

	RTNode* leaf = choose_leaf(stack, entry_idx, stack_size, e, dest_level);
	
	// check if there is space for the new entry
	if (leaf->entry_num < max_entry_num) {
		leaf->entries[leaf->entry_num] = e;//pointer
		leaf->entry_num++;
		/*if (stack_size != 0)
		{
			this->print_node(stack[0],4);
		}	*/
		adjust_tree(stack, entry_idx, stack_size);
		/*if (stack_size != 0)
		{
			this->print_node(stack[0],4);
		}	*/
		delete[] stack;
		delete[] entry_idx;
		return true;
	}

	
	// split is needed.
	bool split = true;
	RTNode* node = leaf;
	Entry new_entry = e;
	while (split) {
		Entry* entry_buffer = new Entry[node->entry_num + 1];
		for (int i = 0; i < node->entry_num; i++) {
			entry_buffer[i] = node->entries[i];
		}
		entry_buffer[max_entry_num] = new_entry;

		int m1, m2;
		linear_pick_seeds(entry_buffer, max_entry_num+1, m1, m2);

		RTNode* new_node = new RTNode(node->level, max_entry_num);
		node->entries[0] = entry_buffer[m1];
		node->entry_num=1;
		new_node->entries[0] = entry_buffer[m2];
		new_node->entry_num = 1;
		// move the selected nodes to the end of the buffer
		swap_entry(entry_buffer, m2, max_entry_num);
		if (m1 == max_entry_num) {
			m1 = m2;
		}
		swap_entry(entry_buffer, m1, max_entry_num-1);
		// (bubble) sort the entries in the remaining set
		// last one should be picked first.
		int remain = max_entry_num-1;
		for (int i = 1; i < remain; i++) {
			for (int j = 0; j < remain - i; j++) {
				if (tie_breaking(entry_buffer[j].get_mbr(), entry_buffer[j+1].get_mbr())) {
					swap_entry(entry_buffer, j, j+1);
				}
			}
		}
		// split procedure
		int max_split_size = (max_entry_num) / 2 + 1;
		BoundingBox old_mbr = node->entries[0].get_mbr();
		BoundingBox new_mbr = new_node->entries[0].get_mbr();
		while (node->entry_num < max_split_size && new_node->entry_num < max_split_size) {
			int old_inc = area_inc(old_mbr, entry_buffer[remain-1].get_mbr());
			int new_inc = area_inc(new_mbr, entry_buffer[remain-1].get_mbr());
			bool add_to_old = false;
			if (old_inc != new_inc) // less enlargement better.
				add_to_old = old_inc < new_inc;
			else if (area(old_mbr) != area(new_mbr)) // smaller area better.
				add_to_old = area(old_mbr) < area(new_mbr);
			else if (node->entry_num != new_node->entry_num) // fewer entries num better.
				add_to_old = node->entry_num < new_node->entry_num;
			else 
				add_to_old = tie_breaking(old_mbr, new_mbr);

			if (add_to_old) {
				node->entries[node->entry_num] = entry_buffer[remain-1];
				node->entry_num++;
				update_mbr(old_mbr, entry_buffer[remain-1].get_mbr());
			}
			else {
				new_node->entries[new_node->entry_num] = entry_buffer[remain-1];
				new_node->entry_num++;
				update_mbr(new_mbr, entry_buffer[remain-1].get_mbr());
			}
			remain--;
		}
		
		// one node reaches max num nodes, assign the remaining to the other node
		if (node->entry_num == max_split_size) {
			for (int i = remain-1; i >= 0; i--) {
				new_node->entries[new_node->entry_num] = entry_buffer[i];
				new_node->entry_num++;
				update_mbr(new_mbr, entry_buffer[i].get_mbr());
			}
		}
		else {
			for (int i = remain-1; i >= 0; i--) {
				node->entries[node->entry_num] = entry_buffer[i];
				node->entry_num++;
				update_mbr(old_mbr, entry_buffer[i].get_mbr());
			}
		}
		// two nodes now. go to a higher level
		if (stack_size == 0) {
			// root reached.
			RTNode* new_root = new RTNode(node->level+1, max_entry_num);
			new_root->entries[0].set_mbr(old_mbr);
			new_root->entries[0].set_ptr(node);
			new_root->entries[1].set_mbr(new_mbr);
			new_root->entries[1].set_ptr(new_node);
			new_root->entry_num = 2;
			root = new_root;
			split = false;
		}
		else {
			stack_size--;
			RTNode* parent = stack[stack_size];
			int idx = entry_idx[stack_size];
			parent->entries[idx].set_mbr(old_mbr);
			new_entry.set_mbr(new_mbr);
			new_entry.set_ptr(new_node);
			if (parent->entry_num < max_entry_num) {
				parent->entries[parent->entry_num] = new_entry;
				parent->entry_num++;
				split = false;
			}
			else
				node = parent;
		}

		delete []entry_buffer;
	}
	adjust_tree(stack, entry_idx, stack_size);

	delete []stack;
	delete []entry_idx;
	return true;
}

bool RTree::del(const vector<int>& coordinate)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}
	/*
	Add your code here
	*/
}


void RTree::query_range(const BoundingBox& mbr, int& result_count, int& node_travelled)
{
	
	result_count = 0;
	node_travelled = 0;
	query_range(root, mbr, result_count, node_travelled);
}


bool RTree::query_point(const vector<int>& coordinate, Entry& result)
{
	BoundingBox mbr(coordinate, coordinate);
	return query_point(root, mbr, result);
}


/**********************************
 *
 * Please do not modify the codes below
 *
 **********************************/

/*********************************************************
  Return true means choose box1 for tie breaking.
  If the two boxes is the same, return true.
  This is to give a unified way of tie-breaking such that if your program is correct, then the result should be same, not influnced by any ties.
 *********************************************************/
bool RTree::tie_breaking(const BoundingBox& box1, const BoundingBox& box2)
{
	//for every dimension, try to break tie by the lowest value, then the highest
	for (int i = 0; i < box1.get_dim(); i++)
	{
		if (box1.get_lowestValue_at(i) != box2.get_lowestValue_at(i))
		{
			return box1.get_lowestValue_at(i) < box2.get_lowestValue_at(i);
		}
		else if (box1.get_highestValue_at(i) != box2.get_highestValue_at(i))
		{
			return box1.get_highestValue_at(i) > box2.get_highestValue_at(i);
		}
	}
	return true;
}


void RTree::stat(RTNode* node, int& record_cnt, int& node_cnt)
{
	if (node->level == 0) {
		record_cnt += node->entry_num;
		node_cnt++;
	}
	else {
		node_cnt++;
		for (int i = 0; i < node->entry_num; i++)
			stat((node->entries[i]).get_ptr(), record_cnt, node_cnt);
	}
}

void RTree::stat()
{
	int record_cnt = 0, node_cnt = 0;
	stat(root, record_cnt, node_cnt);
	cout << "Height of R-tree: " << root->level + 1 << endl;
	cout << "Number of nodes: " << node_cnt << endl;
	cout << "Number of records: " << record_cnt << endl;
	cout << "Dimension: " << dimension << endl;
}


void RTree::print_node(RTNode* node, int indent_level)
{
	BoundingBox mbr = get_mbr(node->entries, node->entry_num);

	char* indent = new char[4*indent_level+1];
	memset(indent, ' ', sizeof(char) * 4 * indent_level);
	indent[4*indent_level] = '\0';

	if (node->level == 0) {
		cout << indent << "Leaf node (level = " << node->level << ") mbr: (";
		for (int i = 0; i < mbr.get_dim(); i++)
		{
			cout << mbr.get_lowestValue_at(i) << " " << mbr.get_highestValue_at(i);
			if (i != mbr.get_dim() - 1)
			{
				cout << " ";
			}
		}
		cout << ")\n";
	}
	else {

		cout << indent << "Non leaf node (level = " << node->level << ") mbr: (";
		for (int i = 0; i < mbr.get_dim(); i++)
		{
			cout << mbr.get_lowestValue_at(i) << " " << mbr.get_highestValue_at(i);
			if (i != mbr.get_dim() - 1)
			{
				cout << " ";
			}
		}
		cout << ")\n";
	}

	Entry *copy = new Entry[node->entry_num];
	for (int i = 0; i < node->entry_num; i++) {
		copy[i] = node->entries[i];
	}

	for (int i = 0; i < node->entry_num; i++) {
		int index = 0; // pick next.
		for (int j = 1; j < node->entry_num - i; j++) {
			if (tie_breaking(copy[j].get_mbr(), copy[index].get_mbr())) {
				index = j;
			}
		}

		if (node->level == 0) {
			Entry& e = copy[index];
			cout << indent << "    Entry: <";
			for (int i = 0; i < e.get_mbr().get_dim(); i++)
			{
				cout << e.get_mbr().get_lowestValue_at(i) << ", ";
			}
			cout << e.get_rid() << ">\n";
		}
		else {
			print_node(copy[index].get_ptr(), indent_level+1);
		}
		// Move the output one to the rear.
		Entry tmp = copy[node->entry_num - i - 1];
		copy[node->entry_num - i - 1] = copy[index];
		copy[index] = tmp;

	}

	delete []indent;
	delete []copy;
}

void RTree::print_tree()
{
	if (root->entry_num == 0)
		cout << "The tree is empty now." << endl;
	else
		print_node(root, 0);
}
