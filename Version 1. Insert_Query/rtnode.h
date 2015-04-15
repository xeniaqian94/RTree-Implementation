#include "boundingbox.h"


class RTNode;

class Entry {
private:
	BoundingBox mbr;
	RTNode* ptr;		//point to the node this entry represents, valid only if this is a non-leaf node entry.
	int rid;			// valid only if this is a leaf node entry.
	
public:
	Entry();
	Entry(const BoundingBox& thatMBR, const int rid);
    Entry(const BoundingBox& thatMBR, RTNode* ptr,const int rid);
	~Entry();
	const BoundingBox& get_mbr() const;
	RTNode* get_ptr() const;
	int get_rid() const;

	void set_mbr(const BoundingBox& thatMBR);
	void set_ptr(RTNode* ptr);
    void set_rid(int rid);

	void print();
};

class RTNode {
	public:
		RTNode(int lev, int size);
        RTNode(int lev, int size, RTNode* parent);
		RTNode(const RTNode& other);
		RTNode& operator=(const RTNode& other);
		~RTNode();
    
        Entry* get_least_enlargement(const BoundingBox& point_box);
        void install(const BoundingBox& point_box, RTNode* ptr, int rid);

	public:
		int entry_num;
		Entry* entries;
		int level;
		int size;
        RTNode* parent;
};
