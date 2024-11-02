#include "Interpreter.h"
#include <stdlib.h>
using namespace std;

void Node::PrintNodeInfo() const
{
	{
		cout << "Node at coordinates [" << coords[1] << "," << coords[0] << "] has a value of: " << value << endl;
		for (Neighbor* checkedNeighbor : neighbors) 
		{
			cout << "It has a neighbor at: [" << checkedNeighbor->neighborNode->coords[1] << "," << checkedNeighbor->neighborNode->coords[0] << "]" << endl;
		}
		cout << endl;
	}
}

bool HashiBoard::ParsePuzzle() 
{
	int rowIndex = 0;
	int columnIndex = 0;
	int currentID = 0;
	for (const vector<int> row : puzzle) 
	{
		for (const int& nodeValue : row)
		{
			if (nodeValue > 0 && nodeValue <= 8)
			{
				int coordinates[2] = { rowIndex, columnIndex };
				islands.push_back(new Node(currentID, nodeValue, coordinates));
				currentID++;
			}
			else if (nodeValue > 8)
			{
				//print warning message that a node greater than 8 is detected
			}
			columnIndex++;
		}
		columnIndex = 0;
		rowIndex++;
	}
	//After parsing, find the neighbors for each island, and print out their info!
	for (Node* checkedNode : islands) 
	{
		UpdateNeighborInfo(checkedNode);
		checkedNode->PrintNodeInfo();
	}
	return true;
}

void HashiBoard::UpdateNeighborInfo(Node* node, bool bShouldClearNeighbors)
{
	Node* neighborNode = nullptr;
	//We need to traverse in each direction in the board in order to find a hit.
	Direction directionsToCheck[4] = { Direction::LEFT, Direction::RIGHT, Direction::DOWN, Direction::UP };
	for (Direction direction : directionsToCheck)
	{
		neighborNode = GetNodeInDirection(direction, node->coords[0], node->coords[1]);
		if (neighborNode)
		{
			//cout << "Node found!" << endl;
			node->neighbors.push_back(new Neighbor(0, direction, neighborNode));
		}
	}
}

Node* HashiBoard::GetNodeAtCoords(int row, int col)
{
	for (Node* checkedNode : islands)
	{
		if (checkedNode->coords[0] == row && checkedNode->coords[1] == col) 
		{
			//cout << "Node found at coords." << endl;
			return checkedNode;
		}
	}
	return nullptr;
}

Node* HashiBoard::GetNodeInDirection(Direction direction, int row, int col)
{
	//Negative values are planned to be for bridges on the graph; check the top of interpreter.h for the commands. 
	int index = 0;
	int parsedNodeValue = 0;
	switch (direction) 
	{
	case Direction::LEFT: 
		for (index = col - 1; index > 0; index--) 
		{
			parsedNodeValue = puzzle[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::RIGHT:
		for (index = col + 1; index < puzzleSizeX; index++)
		{
			parsedNodeValue = puzzle[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::DOWN:
		for (index = row + 1; index < puzzleSizeY; index++)
		{
			parsedNodeValue = puzzle[index][col];
			if (parsedNodeValue > 0)  { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::UP:
		for (index = row - 1; index > 0; index--)
		{
			parsedNodeValue = puzzle[index][col];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	}
	return nullptr;
}