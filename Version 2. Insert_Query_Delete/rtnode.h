#include "boundingbox.h"
#include <vector>


class RTNode;

class Entry {
private:
	BoundingBox mbr;
	RTNode* ptr;		//point to the node this entry represents, valid only if this is a non-leaf node entry.
	int rid;			// valid only if this is a leaf node entry.
	
public:
	Entry();
	Entry(const BoundingBox& thatMBR, const int rid);
	~Entry();
	const BoundingBox& get_mbr() const;
	RTNode* get_ptr() const;
	int get_rid() const;

	void set_mbr(const BoundingBox& thatMBR);
	void set_ptr(RTNode* ptr);

	void print();
};

class RTNode {
	public:
		RTNode(int lev, int size);        
		RTNode(const RTNode& other);
		RTNode& operator=(const RTNode& other);
		~RTNode();

	public:
		int entry_num;
		Entry* entries;
		int level;
		int size;
};
