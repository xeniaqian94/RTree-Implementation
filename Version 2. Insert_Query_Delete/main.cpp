#include <cstdlib>
#include <fstream>
#include "rtree.h"

using namespace std;

const int MAX_CMD_LEN = 256;
const int DOMAIN_SIZE = 10000;

void help()
{
	cout << "============================================================================\n";
	cout << "Commands:\n";
	cout << "============================================================================\n";
	cout << "i x1(int) x2(int) ... xd(int) rid(int) : insert a record with d-dimension key (x1, x2,... , xd) and record id rid\n";
	cout << "d x1(int) x2(int) ... xd(int) : delete the record with key (x1, x2,... , xd)\n";
	cout << "ri s(int) num(int) : random insertions of num records with seed s\n";
	cout << "rd s(int) num(int) : random deletions of num records with seed s\n";
	cout << "qp x1(int) x2(int) ... xd(int) : query the record with key (x1, x2, ... , xd)\n";
	cout << "qr x1min(int) x1max(int) x2min(int) x2max(int) ... xdmin(int) xdmax(int) : find records inside range\n";
	cout << "     where ximin<=xi<=ximax\n";
	cout << "s : print the statistic information of the tree\n";
	cout << "p : print the tree\n";
	cout << "h : show this help menu\n";
	cout << "x : exit\n";
	cout << "============================================================================\n";
}

void error(const char* cmd)
{
	cerr << "Error: " << cmd << endl;
}

bool process(char* cmd, RTree& tree, int dimension)
{
	
	const int MAX_ARG_NUM = 256; // limit to at most 256 arguments
	char* args[MAX_ARG_NUM];
	
	char msg[1024]; // error message.
	int actualMaxArgNum = 1 + dimension * 2;
	if (actualMaxArgNum > MAX_ARG_NUM)
	{
		sprintf(msg, "Too many command arguments");
		error(msg);
		return true;
	}
	int num_arg = 0;
	char* token = strtok(cmd, "\n \t");
	for (; token != NULL && num_arg < actualMaxArgNum; num_arg++) {
		args[num_arg] = token;
		token = strtok(NULL, " \t");
	}
	if (num_arg == 0 || token != NULL) {
		sprintf(msg, "Wrong number of command arguments");
		error(msg);
		return true;
	}
	if (strcmp(args[0], "i") == 0) { // insertion.
		if (num_arg != dimension + 2) {
			sprintf(msg, "Wrong number of arguments for command 'i'");
			error(msg);
		}
		else {
			//insert a point, modelled by a bounding box
			vector<int> coordinate;
			for (int i = 0; i < dimension; i++)
			{
				int coord = atoi(args[i + 1]);
				coordinate.push_back(coord);
			}
			int rid = atoi(args[dimension + 1]);
			try {
				if (tree.insert(coordinate, rid))
					cout << "Insertion done.\n";
				else
					cout << "Insertion failed.\n";
			}
			catch (bad_alloc& ba)  {
				sprintf(msg, "bad_alloc caught <%s> ", ba.what());
				error(msg);
			}
		}
		return true;
	}
	else if (strcmp(args[0], "d") == 0) { // deletion.
		if (num_arg != dimension + 1) {
			sprintf(msg, "Wrong number of arguments for command 'd'");
			error(msg);
		}
		else {
			vector<int> coordinate;
			for (int i = 0; i < dimension; i++)
			{
				int coord = atoi(args[i + 1]);
				coordinate.push_back(coord);
			}

			if (tree.del(coordinate))
				cout << "Deletion done.\n";
			else
				cout << "Deletion failed.\n";
		}
		return true;
	}
	else if (strcmp(args[0], "ri") == 0) { // random insertion.
		if (num_arg != 3) {
			sprintf(msg, "Wrong number of arguments for command 'ri'");
			error(msg);
		}
		else {
			srand(atoi(args[1]));
			int num = atoi(args[2]);
			int succeed = 0;
			for (int i = 0; i < num; i++) {
				vector<int> coordinate;
				for (int j = 0; j < dimension; j++)
				{
					int coord = rand() % DOMAIN_SIZE;
					coordinate.push_back(coord);
				}
				int rid = rand();
				
				try {
					if (tree.insert(coordinate, rid)) {
						succeed++;
					}
				}
				catch (bad_alloc& ba)  {
					sprintf(msg, "bad_alloc caught <%s> ", ba.what());
					error(msg);
				}
				//tree.print_tree();
			}
			cout << succeed << " out of " << num << " insertion(s) suceeded.\n";
		}
		return true;
	}
	else if (strcmp(args[0], "rd") == 0) { // random deletion.
		if (num_arg != 3) {
			sprintf(msg, "Wrong number of arguments for command 'rd'");
			error(msg);
		}
		else {
			srand(atoi(args[1]));
			int num = atoi(args[2]);
			int succeed = 0;
			for (int i = 0; i < num; i++) {
				vector<int> coordinate;
				for (int j = 0; j < dimension; j++)
				{
					int coord = rand() % DOMAIN_SIZE;
					coordinate.push_back(coord);
				}
				int dummy = rand(); // to be compatible with ``ri''.
				if (tree.del(coordinate)) {
					succeed++;
				}
			}
			cout << succeed << " out of " << num << " deletion(s) suceeded.\n";
		}
		return true;
	}
	else if (strcmp(args[0], "qr") == 0) { // range query.
		if (num_arg != 1 + dimension * 2) {
			sprintf(msg, "Wrong number of arguments for command 'qr'");
			error(msg);
		}
		else {
			vector<int> lowest;
			vector<int> highest;
			for (int i = 0; i < dimension; i++)
			{
				lowest.push_back(atoi(args[1 + i*2]));			
				highest.push_back(atoi(args[2 + i*2]));
			}

			BoundingBox mbr(lowest, highest);

			int result_count = 0;
			int node_travelled = 0;
			tree.query_range(mbr, result_count, node_travelled);
			cout << "Number of results: " << result_count << endl;
			cout << "Number of nodes visited: " << node_travelled << endl;
		}
		return true;
	}
	else if (strcmp(args[0], "qp") == 0) { // point query.
		if (num_arg != 1 + dimension) {
			sprintf(msg, "Wrong number of arguments for command 'qp'");
			error(msg);
		}
		else {
			Entry result;

			vector<int> coordinate;

			for (int i = 0; i < dimension; i++)
			{
				coordinate.push_back(atoi(args[i + 1]));
			}

			if (tree.query_point(coordinate, result)) {
				cout << "Record: <";
				BoundingBox resultP = result.get_mbr();
				for (int i = 0; i < resultP.get_dim(); i++)
				{
					cout << resultP.get_lowestValue_at(i);
					if (i != resultP.get_dim() - 1)
					{
						cout << ", ";
					}
				}
				cout << ", " << result.get_rid()  << ">\n";
			}
			else {
				cout << "Record not found.\n";
			}
		}
		return true;
	}
	else if (strcmp(args[0], "s") == 0) { // statistics.
		tree.stat();
		return true;
	}
	else if (strcmp(args[0], "p") == 0) { // print tree.
		tree.print_tree();
		return true;
	}
	else if (strcmp(args[0], "h") == 0) { // print help menu.
		help();
		return true;
	}
	else if (strcmp(args[0], "x") == 0) { // exit
		return false;
	}
	else {
		sprintf(msg, "Invalid command '%s'.\nType 'h' to print the help menu.", args[0]);
		error(msg);
		return true;
	}
}



int main(int argc, char *argv[])
{//argc also counts the argv[0] that is the name of the program
	
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " Max_#entries_in_a_node Dimensionality_of_Rtree\" [file_containing_commmands].\n";
		return 0;
	}

	// Create an R-tree.
	int max_entry_num = atoi(argv[1]);
	if (max_entry_num < 2) {
		cerr << "Number of entries should be an integer > 2.\n";
		return 0;
	}
	int dimension = atoi(argv[2]);
	RTree tree(max_entry_num, dimension);

	// Processing input commands.
	char command[MAX_CMD_LEN];
	if (argc == 4) {
		ifstream fin(argv[3]);
		while (fin.getline(command, MAX_CMD_LEN)) {
			//cout << command << endl;
			if (! process(command, tree, dimension))
				break;
		}
	}
	else {
		while (true) {
			cout << ">> ";
			cin.getline(command, MAX_CMD_LEN);
			if (! process(command, tree, dimension))
				break;
		}
	}

	return 0;
}
