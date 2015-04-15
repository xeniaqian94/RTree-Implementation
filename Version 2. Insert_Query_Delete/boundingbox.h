
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;


class BoundingBox {
private:
	vector<int> lowest; //lowest coordinate of the bounding box
	vector<int> highest; //highest coordinate
public:
	BoundingBox();
	BoundingBox(vector<int> thatLow, vector<int> thatHigh);
	BoundingBox(const BoundingBox& thatBox);

	const vector<int>& get_lowest() const;
	const vector<int>& get_highest() const;
	int get_dim() const;
	int get_area() const;
	int get_lowestValue_at(const int index) const;
	int get_highestValue_at(const int index) const;

	bool is_equal(const BoundingBox& rhs) const; // if this mbr equals to rhs mbr
	bool is_intersected(const BoundingBox& rhs) const;// if this mbr overlaps with rhs mbr
	bool is_valid() const;
	void print() const;

	void group_with(const BoundingBox& rhs); //update this by the MBR of this and rhs
	void set_boundingbox(const BoundingBox& rhs);
};

