/* Implementations of R tree */
#include <cmath>
#include "rtree.h"


const double EPSILON = 1E-10;
const int DOMAIN_SIZE = 10000;

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


bool RTree::insert(const vector<int>& coordinate, int rid)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}
    Entry check_entry;
    if (query_point(coordinate, check_entry)) return false;
    
    BoundingBox point_box(coordinate,coordinate);
    RTNode* result_leaf=choose_leaf(point_box, rid);
    RTNode* result_leaf_ll=NULL;
    if (result_leaf->entry_num<result_leaf->size){ // no SPLIT! but need to adjust all the levels of mbr
        result_leaf->install(point_box,NULL,rid);
    }
    
    else {
//        cout<<"leaf is full causing SPLIT!"<<endl;
        Entry* entry_list = new Entry[result_leaf->entry_num+1];
        
        for (int i = 0; i < result_leaf->entry_num; i++)
            entry_list[i] = result_leaf->entries[i];
        Entry e(point_box,rid);
        entry_list[result_leaf->entry_num]=e;
        RTNode* old_full_leaf=result_leaf;
        split_node(result_leaf,result_leaf_ll,entry_list); //split M+1 entries into 2 nodes
        if (old_full_leaf->parent!=NULL){
            RTNode* parent=old_full_leaf->parent;
            for (int i=0;i<parent->entry_num;i++)
                if (parent->entries[i].get_ptr()==old_full_leaf){
                    parent->entries[i].set_ptr(result_leaf); delete old_full_leaf;break;}
        }
    }
    RTNode* p;
    RTNode* pp=NULL;
    adjust_tree(result_leaf,result_leaf_ll,p,pp);
    if (pp!=NULL){       //if the root was split
        RTNode* new_root=new RTNode(p->level+1,p->size,NULL);
        new_root->install(get_mbr(p->entries, p->entry_num),p,-1);
        new_root->install(get_mbr(pp->entries,pp->entry_num),pp,-1);
        root=new_root;
    }
    else root=p;
    return true;
}


RTNode* RTree::choose_leaf(const BoundingBox& point_box, int rid)
{
    
    RTNode* N=root;
    if (N->entry_num==0) return root;
    else while(N->entries[0].get_ptr()!=NULL){ 
        Entry* F=N->get_least_enlargement(point_box);//find Entry F whose mbr needs least enlargement to include E
        N=F->get_ptr();
        
    }
    return N; //a leaf node
}

void RTree::adjust_tree(RTNode*& L, RTNode*& LL,RTNode*& p,RTNode*& pp)
{
    RTNode* N=L;
    RTNode* NN=NULL;
    if (LL!=NULL) NN=LL;
    RTNode* old_full_p=NULL;
    while (N->parent!=NULL)
    {
        p=N->parent;
        for (int i=0;i<p->entry_num;i++){
            if (p->entries[i].get_ptr()==N){ //no split on the last recursion
                p->entries[i].set_mbr(get_mbr(N->entries,N->entry_num));
                break;
            }
            if (old_full_p!=NULL){ //split happened on the last recursion
                if (p->entries[i].get_ptr()==old_full_p){
                    p->entries[i].set_ptr(N);
                    p->entries[i].set_mbr(get_mbr(N->entries,N->entry_num));
                    break;
                        
                }
            }
        }
        if (NN!=NULL){
            if (p->entry_num<p->size){
                NN->parent=p;
                p->install(get_mbr(NN->entries,NN->entry_num),NN,-1);
                pp=NULL;
                old_full_p=NULL;
            }
            else if (p->entry_num==p->size){ 
                Entry* entry_list = new Entry[p->entry_num+1];
                for (int i = 0; i < p->entry_num; i++)
                    entry_list[i] = p->entries[i];
                Entry ENN(get_mbr(NN->entries,NN->entry_num),NN,-1);
                entry_list[p->entry_num]=ENN;
                old_full_p=p;
                split_node(p,pp,entry_list);
            }
        }
        N=p; NN=pp;
    }
    p=N;pp=NN;
}


void RTree::split_node(RTNode*& leaf_l, RTNode*& leaf_ll,Entry* entry_list)
{
    int seed_x=0, seed_y=0;
    linear_pick_seed(entry_list, seed_x, seed_y);
    
    RTNode* L=new RTNode(leaf_l->level,leaf_l->size,leaf_l->parent); //inherits the original parent
    RTNode* LL=new RTNode(leaf_l->level,leaf_l->size,leaf_l->parent); //no default parent
    L->install(entry_list[seed_x].get_mbr(),entry_list[seed_x].get_ptr(),entry_list[seed_x].get_rid());
    
    BoundingBox bL(entry_list[seed_x].get_mbr());
    
    LL->install(entry_list[seed_y].get_mbr(),entry_list[seed_y].get_ptr(),entry_list[seed_y].get_rid());
    BoundingBox bLL(entry_list[seed_y].get_mbr());
    
    int* flag=new int[max_entry_num+1];
    flag[seed_x]=1;
    flag[seed_y]=1;
    
    while (L->entry_num<(int)((max_entry_num+2)/2)&&LL->entry_num<(int)((max_entry_num+2)/2))
    {
        int next_index=-1;
        int group=1;
        for (int i=0;i<max_entry_num+1;i++)
        {
            if (flag[i]!=1){ //not picked yet
                if (next_index==-1) next_index=i;
                else if (tie_breaking(entry_list[i].get_mbr(), entry_list[next_index].get_mbr()))
                    next_index=i;
            }
        }
        BoundingBox temp_L(bL);
        temp_L.group_with(entry_list[next_index].get_mbr());
        
        BoundingBox temp_LL(bLL);
        temp_LL.group_with(entry_list[next_index].get_mbr());
        
        int d1=temp_L.get_area()-bL.get_area();
        
        int d2=temp_LL.get_area()-bLL.get_area();
        
        if ((d1-d2)>0) group=2;
        else if (d1==d2){
            if (bL.get_area()<bLL.get_area()) group=1;
            else if (bL.get_area()>bLL.get_area()) group=2;
            else{
                if (L->entry_num<LL->entry_num) group=1;
                else if (L->entry_num>LL->entry_num) group=2;
                else{
                    if (tie_breaking(bL,bLL)) group=1;
                    else group=2;
                }
             
            }
        }
        else if ((d1-d2)<0) group=1;
        
        flag[next_index]=1;
        if (group==1){
            L->install(entry_list[next_index].get_mbr(),entry_list[next_index].get_ptr(),entry_list[next_index].get_rid());
            bL.group_with(entry_list[next_index].get_mbr());
        }
        else{
            LL->install(entry_list[next_index].get_mbr(),entry_list[next_index].get_ptr(),entry_list[next_index].get_rid());
            bLL.group_with(entry_list[next_index].get_mbr());
        }
    }
    if (L->entry_num<LL->entry_num){
        for (int i=0;i<max_entry_num+1;i++)
        {
            if (i!=seed_x&&i!=seed_y&&flag[i]!=1){ //pick_next&install it
                L->install(entry_list[i].get_mbr(),entry_list[i].get_ptr(),entry_list[i].get_rid());
            }
        }
    }
    else if (LL->entry_num<L->entry_num){
        for (int i=0;i<max_entry_num+1;i++)
        {
            if (i!=seed_x&&i!=seed_y&&flag[i]!=1){ //pick_next&install it
                LL->install(entry_list[i].get_mbr(),entry_list[i].get_ptr(),entry_list[i].get_rid());
            }
        }
    }
    leaf_l=L;
    leaf_ll=LL;
}

void RTree::linear_pick_seed(Entry* entry_list,int& seed_x, int& seed_y)
{
    double max_normalized_separation=0;
    for (int i=0;i<dimension;i++){
        find_pairs(i,entry_list,max_normalized_separation,seed_x,seed_y);
    }
    if (seed_x==seed_y){
        int temp=0; //the smallest entry index after sort
        if (seed_x==0) temp=1;
        for (int i=0;i<max_entry_num+1;i++)
            if (i!=seed_x && tie_breaking(entry_list[temp].get_mbr(),entry_list[i].get_mbr())) temp=i;
        seed_y=temp;
    }
    
}

void RTree::find_pairs(int dimension, Entry* entry_list, double& max, int& x, int &y)
{
    int highest_low_side=0,lowest_high_side=DOMAIN_SIZE;
    int lowest_low_side=DOMAIN_SIZE,highest_high_side=0;
    int this_y=0,this_x=0;
    for (int i=0;i<max_entry_num+1;i++){
        BoundingBox this_box(entry_list[i].get_mbr());
        int this_low=this_box.get_lowestValue_at(dimension);
        int this_high=this_box.get_highestValue_at(dimension);
        if (this_low>highest_low_side){
            highest_low_side=this_low;
            this_y=i;
        }
        else if (this_low==highest_low_side&&tie_breaking(this_box,entry_list[this_y].get_mbr()))
            this_y=i;
    
        if (this_high<lowest_high_side){
            lowest_high_side=this_high;
            this_x=i;
        }
        else if (this_high==lowest_high_side&&tie_breaking(this_box, entry_list[this_x].get_mbr()))
            this_x=i;
        if (this_low<lowest_low_side) lowest_low_side=this_low;
        if (this_high>highest_high_side) highest_high_side=this_high;
    }
    
    double normalized_separation=1.00*(highest_low_side-lowest_high_side)/(highest_high_side-lowest_low_side);
//    cout<<this_y<<","<<this_x<<"normalized_separation"<<normalized_separation<<endl;
    if (normalized_separation>max){
        y=this_y;
        x=this_x;
        max=normalized_separation;
    }
    
}

void RTree::query_range(const BoundingBox& mbr, int& result_count, int& node_travelled)
{
	if (mbr.get_dim() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}
    if (root->entry_num==0) node_travelled=1;
    else query_range_subtree(root,mbr,result_count,node_travelled);
}
void RTree::query_range_subtree(RTNode* root, const BoundingBox& mbr, int& result_count, int& node_travelled)
{
    RTNode* N=root;
    node_travelled++;
    if (N->entries[0].get_ptr()!=NULL){ //not a leaf
        for (int i=0;i<N->entry_num;i++)
            if (N->entries[i].get_mbr().is_intersected(mbr))
                query_range_subtree(N->entries[i].get_ptr(), mbr, result_count,node_travelled);
    }
    else if (N->entries[0].get_ptr()==NULL){
        for (int i=0;i<N->entry_num;i++)
            if (N->entries[i].get_mbr().is_intersected(mbr))
                result_count++;
    }
}


bool RTree::query_point(const vector<int>& coordinate, Entry& result)
{
	if (coordinate.size() != this->dimension)
	{
		cerr << "R-tree dimensionality inconsistency\n";
	}
    if (root->entry_num==0) return false;
    else return query_point_subtree(root,coordinate,result);
}
bool RTree::query_point_subtree(RTNode* root, const vector<int>& coordinate, Entry& result)
{
    RTNode* N=root;
    BoundingBox b(coordinate,coordinate);
    if (N->entries[0].get_ptr()!=NULL){ //not a leaf
        for (int i=0;i<N->entry_num;i++){
            
            if (N->entries[i].get_mbr().is_intersected(b))
                if (query_point_subtree(N->entries[i].get_ptr(), coordinate, result)) return true;
        }
    }
    else if (N->entries[0].get_ptr()==NULL){ //is a leaf
        for (int i=0;i<N->entry_num;i++){
//            N->entries[i].get_mbr().print();
            if (N->entries[i].get_mbr().is_equal(b)){
                result=N->entries[i];
                return true;
            }
        }

    }
    
    return false;
}


/**********************************
 *
 * Please do not modify the codes below
 *
 **********************************/

//
// Calcuate the MBR of a set of entries, of size ``len''.
// Store the MBR in the first entry
//
BoundingBox RTree::get_mbr(Entry* entry_list, int len)
{
	BoundingBox mbr(entry_list[0].get_mbr());
	for (int i = 1; i < len; i++) {        
		mbr.group_with(entry_list[i].get_mbr());
	}
	return mbr;
}


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

	if (node->level == 0) {   //leaf node
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
