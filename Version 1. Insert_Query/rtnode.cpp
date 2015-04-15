#include "rtnode.h"

//======================== Entry implementation =====================================================
const int DOMAIN_SIZE = 10000;

Entry::Entry():mbr() {
	this->rid = -1;
	this->ptr = NULL;
}

Entry::Entry(const BoundingBox& thatMBR, const int rid):mbr(thatMBR) {
    this->rid = rid;
	this->ptr = NULL;
}

Entry::Entry(const BoundingBox& thatMBR, RTNode* ptr,const int rid):mbr(thatMBR){
    this->rid = rid;
    this->ptr = ptr;
}

Entry::~Entry() {
	this->ptr = NULL;
}

const BoundingBox& Entry::get_mbr() const {
	return this->mbr;
}	


RTNode* Entry::get_ptr() const {
	return this->ptr;
}

int Entry::get_rid() const { 
	return this->rid;
}


void Entry::set_mbr(const BoundingBox& thatMBR) {
	this->mbr.set_boundingbox(thatMBR);//ATTENTION: should be deep copy...
}

void Entry::set_rid(int rid) {
    this->rid=rid;
}

void Entry::set_ptr(RTNode* ptr) {
	this->ptr = ptr;
}


void Entry::print() {
	this->mbr.print();
	cout << this->rid << endl;
	cout << this->ptr << endl;
}

//======================== RTNode implementation ==============================================

RTNode::RTNode(int lev, int s)
{
	entry_num = 0;
	entries = new Entry[s];
	level = lev;
	size = s;
    this->parent=NULL;
}
RTNode::RTNode(int lev, int s,RTNode* parent)
{
    entry_num = 0;
    entries = new Entry[s];
    level = lev;
    size = s;
    this->parent=parent;
}

RTNode::RTNode(const RTNode& other)
{
	entries = new Entry[other.size];
	*this = other;
}

RTNode& RTNode::operator=(const RTNode& other)
{
	if (&other != this) {
		entry_num = other.entry_num;
		level = other.level;
		size = other.size;
		for (int i = 0; i < entry_num; i++)
			entries[i] = other.entries[i];
	}
	return *this;
}


RTNode::~RTNode()
{
	if (level != 0) {
		for (int i = 0; i < entry_num; i++) {
            delete entries[i].get_ptr();
			entries[i].set_ptr(NULL);
		}
	}
	delete []entries;
	entries = NULL;
}

Entry* RTNode::get_least_enlargement(const BoundingBox& point_box){
    int least_index=0;
    int least_enlargement_area=DOMAIN_SIZE*DOMAIN_SIZE;
    int least_rectangle_area=DOMAIN_SIZE*DOMAIN_SIZE;
    for (int i=0;i<entry_num;i++){
        BoundingBox temp(entries[i].get_mbr());
        int this_area=temp.get_area();
        temp.group_with(point_box);
        int new_area=temp.get_area();
        int enlargement_area=new_area-this_area;
        if (enlargement_area<least_enlargement_area) {
            least_index=i;
            least_enlargement_area=enlargement_area;
            least_rectangle_area=this_area;
        }
        else if (enlargement_area==least_enlargement_area&&this_area<least_rectangle_area) {
            least_index=i;
            least_rectangle_area=this_area;
        }
        else if (enlargement_area==least_rectangle_area&&this_area==least_rectangle_area){
            if (temp.tie_breaking(entries[i].get_mbr(),entries[least_index].get_mbr())){
                least_index=i;
                least_rectangle_area=this_area;
            }
        }
    }
    return (entries+least_index);
}

void RTNode::install(const BoundingBox& point_box, RTNode* ptr, int rid)
{
    if (entry_num<size){
        Entry e(point_box,ptr,rid);
        entries[entry_num]=e;
        if (ptr!=NULL) ptr->parent=this;
        entry_num++;
    }
    
        
}


