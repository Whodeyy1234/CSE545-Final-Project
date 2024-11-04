#include <vector>
#include <iostream>
#include <SDL2Singleton.h>
#pragma once
using namespace std;

struct Neighbor;

// Helpful macros for rendering.
#define VERT_SINGLE_BRIDGE -1
#define VERT_DOUBLE_BRIDGE -2
#define HORI_SINGLE_BRIDGE -3
#define HORI_DOUBLE_BRIDGE -4
#define GRID_HEIGHT(y) static_cast<int>(SCREEN_HEIGHT / y)
#define GRID_WIDTH(x) static_cast<int>(SCREEN_WIDTH / x)
#define ISLAND_RADIUS_FACTOR 1/6
#define OFFSET_FACTOR 1/4
#define BRIDGE_LENGTH_FACTOR 3/4

enum class Direction 
{
	INVALID,
	LEFT,
	DOWN,
	RIGHT,
	UP,
	MAX
};

//The baseline structure, these are the 'islands' in the graph. 
class Node
{
public:
	int nodeID; //ID of the node for quick comparison purposes
	int value; //The value of the node, that is, how many bridges should the node have?
	int coords[2]; //The coordinates of the node, which are in (y, x) coordinates due to how arrays are structured.
	bool bIsComplete; //Is the node completely filled with bridges? If so, should disregard from connection criteria for other bridges. 

	Node(int argNodeID, int argValue, int argCoords[2])
	{
		nodeID = argNodeID;
		value = argValue;
		coords[0] = argCoords[0];
		coords[1] = argCoords[1];
		bIsComplete = false;
	}
	~Node()
	{
		for (struct Neighbor* neighbor : neighbors)
		{
			delete neighbor;
		}
	}

	void PrintNodeInfo() const;

	//A 'neighbor' is a struct that contains the # of bridges between the two nodes, the direction to the neighbor, as well as a pointer to the neighboring node object itself. 
	//Two nodes who are neighbors should have pointers to one another. 
	vector<Neighbor*> neighbors;
};

/// <summary>
/// A 'neighbor' is a struct that contains the # of bridges between the two nodes, the direction to the neighbor, as well as a pointer to the neighboring node object itself. 
/// Two nodes who are neighbors should have pointers to one another. 
/// </summary>
struct Neighbor
{
	int numOfBridges; //Should ONLY be 0, 1, or 2. 
	Direction neighborDirection; //The direction that the neighbor is located towards. Should help to easily build bridges in the specified direction. 
	Node* neighborNode; //Actual pointer to the node of the neighbor.
	Neighbor(int argBridges, Direction argNeighborDirection, Node* neighborNodePointer)
	{
		numOfBridges = argBridges;
		neighborDirection = argNeighborDirection;
		neighborNode = neighborNodePointer;
	}
};

/// <summary>
/// The actual board itself, only one should be generated from input when the program is run. 
/// </summary>
class HashiBoard
{
public:
	//TO-DO: Add a way to handle setting up a puzzle through a file, probably a txt. 

	const int puzzleSizeX = 7;
	const int puzzleSizeY = 7;

	//A vector of vectors that indicates the state of the puzzle. Each element of the puzzle dictates a different part of the puzzle based on their value.
	//0:	An empty spacee on the puzzle. 
	//1-8:	An island on the puzzle, with a value as dictated with the number.
	//-1:	|	1 bridge vertical
	//-2:	||	2 bridges vertical
	//-3:	-	1 bridge horizontal
	//-4:	=	2 bridges horizontal
	vector<vector<int>> puzzle =
	{
		{0, 2, 0, 5, 0, 0, 2 },
		{0, 0, 0, 0, 0, 0, 0 },
		{4, 0, 2, 0, 2, 0, 4 },
		{0, 0, 0, 0, 0, 0, 0 },
		{0, 1, 0, 5, 0, 2, 0 },
		{0, 0, 0, 0, 0, 0, 0 },
		{4, 0, 0, 0, 0, 0, 3 }
	};

	vector<Node*> islands; //All the islands that are on the board.

	/// <summary>
	/// Default constructor to help with initialization.
	/// </summary>
	HashiBoard() 
	{
		texture = SDL_CreateTexture(SDL->GetRenderer(), SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, SCREEN_WIDTH, SCREEN_HEIGHT);
		bLongerWidth = puzzle[0].size() >= puzzle.size();
	}
	/// <summary>
	/// Default destructor to help free memory.
	/// </summary>
	~HashiBoard()
	{
		SDL_DestroyTexture(texture);
		for (Node* node : islands)
		{
			delete node;
		}
	}

	/// <summary>
	/// Parses through a given puzzle in the puzzle variable, and converts it to islands, and figures out the neighbors to the islands. 
	/// </summary>
	/// <returns>Returns true if puzzle was successfully setup, false otherwise. </returns>
	bool ParsePuzzle();

	/// <summary>
	/// Update a given node's neighbors by checking for the first 'hit' of a node in each cardinal direction.
	/// If a 'hit' gets nothing, or hits a bridge, then no neighbor exists in that direction. 
	/// </summary>
	/// <param name="node">Pointer to the node that is checking its neighbors.</param>
	/// <param name="bShouldClearNeighbors">Should the neighbors list be cleared out before updating? Should usually be true, but here just in case. </param>
	void UpdateNeighborInfo(Node* node, bool bShouldClearNeighbors = true); 

	/// <summary>
	/// Gets the node at specified (y, x) coordinates. Y
	/// Yes, I hate that they're in (y,x), but at least I'm making it clear here.
	/// </summary>
	/// <param name="row"></param>
	/// <param name="col"></param>
	/// <returns>Returns a valid node pointer if one exists, otherwise returns nullptr. </returns>
	Node* GetNodeAtCoords(int row, int col);

	/// <summary>
	/// Checks if a node exists in any direction from a given point, and returns it. 
	/// </summary>
	/// <param name="direction"></param>
	/// <param name="row"></param>
	/// <param name="col"></param>
	/// <returns>Returns a valid node pointer if one exists, otherwise returns nullptr.</returns>
	Node* GetNodeInDirection(Direction direction, int row, int col);

	void RenderBoard();

private:
	/// <summary>
	/// Texture to render to.
	/// </summary>
	SDL_Texture* texture;
	/// <summary>
	/// Whether or not the width is larger than the height.
	/// </summary>
	bool bLongerWidth;

	/// <summary>
	/// Obtains the radius of the islands.
	/// </summary>
	/// <returns>Radius of the islands.</returns>
	const int GetIslandRadius() const;
	/// <summary>
	/// Renders a single island to the texture.
	/// </summary>
	/// <param name="Id">Id of the given island.</param>
	/// <param name="CenterX">The center of the island in pixels x.</param>
	/// <param name="CenterY">The center of the island in pixels y.</param>
	void RenderIsland(const int Id, const int CenterX, const int CenterY);
	/// <summary>
	/// Renders a single bridge connection to the texture.
	/// </summary>
	/// <param name="Type">The type of bridge to build.</param>
	/// <param name="CenterX">The center of the island in pixels x.</param>
	/// <param name="CenterY">The center of the island in pixels y.</param>
	void RenderBridge(const int Type, const int CenterX, const int CenterY);
};