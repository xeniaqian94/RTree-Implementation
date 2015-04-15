#include "boundingbox.h"

//======================== BoundingBox implementation ==============================================
BoundingBox::BoundingBox() {
}

BoundingBox::BoundingBox(vector<int> thatLow, vector<int> thatHigh) {
	if (thatHigh.size() != thatLow.size())
	{
		cerr << "lowest and highest point of rectangle should have the same length\n";
		//exit(-1);
	}

	this->lowest = thatLow;
	this->highest = thatHigh;
	this->is_valid();
}

BoundingBox::BoundingBox(const BoundingBox& thatBox) {
	this->lowest = thatBox.lowest;
	this->highest = thatBox.highest;
}

const vector<int>& BoundingBox::get_lowest() const {
	return this->lowest;
}

const vector<int>& BoundingBox::get_highest() const {
	return this->highest;
}

int BoundingBox::get_dim() const {
	return (int)this->lowest.size();
}

int BoundingBox::get_area() const {
	int area = 1;

	for (int cIndex = 0; cIndex < this->get_dim(); cIndex++)
	{
		area *= this->highest[cIndex] - this->lowest[cIndex];
	}

	return area;
}

int BoundingBox::get_lowestValue_at(const int index) const {
	return this->lowest[index];
}

int BoundingBox::get_highestValue_at(const int index) const {
	return this->highest[index];
}

//if two bounding boxes are the same with respect to their coordinates
bool BoundingBox::is_equal(const BoundingBox& rhs) const {
	if (this->lowest == rhs.get_lowest() && this->highest == rhs.get_highest())
	{
		return true;
	}
	return false;
}

bool BoundingBox::is_intersected(const BoundingBox& rhs) const { 
	
	if (this->get_dim() != rhs.get_dim())
	{
		cerr << "domensiomality inconsistency" << endl;
	}

	const vector<int>& thatLow = rhs.get_lowest();
	const vector<int>& thatHigh = rhs.get_highest();

	//if the two shapes intersect, they must intersect in all dimensions.
	for (int cIndex = 0; cIndex < this->get_dim(); cIndex++)
	{
		if (this->lowest[cIndex] > thatHigh[cIndex] || this->highest[cIndex] < thatLow[cIndex]) return false;
	}
	return true;
}

bool BoundingBox::is_valid() const {
	for (int i = 0; i < this->get_dim(); i++)
	{
		if (this->lowest[i] > this->highest[i])
		{
			cerr << "bounding box has low value " << this->lowest[i] << " larger than high" << this->highest[i] <<" at dimension " << i << ", which should not\n";
			//exit(-1);
			return false;
		}
	}
	return true;
}

void BoundingBox::print() const {
	cout << "bounding box (";
	for (int i = 0; i < this->get_dim(); i++)
	{
		cout << this->lowest[i];
		if (i != this->get_dim() - 1)
		{
			cout << " ";
		}
	}
	cout << ",";
	for (int i = 0; i < this->get_dim(); i++)
	{
		cout << this->highest[i];
		if (i != this->get_dim() - 1)
		{
			cout << " ";
		}
	}
	cout << ")\n";
}


void BoundingBox::group_with(const BoundingBox& rhs) {
	if (this->get_dim() != rhs.get_dim()) {
		cerr << "domensiomality inconsistency" << endl;
		//exit(-1);
	};

	const vector<int>& thatLow = rhs.get_lowest();
	const vector<int>& thatHigh = rhs.get_highest();

	for (int cIndex = 0; cIndex < this->get_dim(); cIndex++)
	{
		this->lowest[cIndex] = this->lowest[cIndex] <= thatLow[cIndex] ? this->lowest[cIndex] : thatLow[cIndex];
		this->highest[cIndex] = this->highest[cIndex] >= thatHigh[cIndex] ? this->highest[cIndex] : thatHigh[cIndex];
	}
}

void BoundingBox::set_boundingbox(const BoundingBox& rhs) {

	this->lowest = rhs.get_lowest();
	this->highest = rhs.get_highest();
}

bool BoundingBox::tie_breaking(const BoundingBox& box1, const BoundingBox& box2)
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

 