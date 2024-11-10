#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <SDL2Singleton.h>
#include <unordered_set>
#include <stack>
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
#define BITMASK_BOUNDARY 4

enum class Direction : int
{
	INVALID = -1,
	LEFT = 0,
	RIGHT = 1,
	DOWN = 2,
	UP = 3,
	MAX
};

/// <summary>
/// Parameters used for the algorithm.
/// </summary>
struct Parameters
{
	unsigned int seed = 0;
	int populationSize = 0;
	float crossoverProb = 0.f;
	float mutationProb = 0.f;
	int maxGenerations = 0;
	bool bWithWisdom = false;
	int gensPerWisdom = 0;
	float elitismPerc = 0.f;
};

struct BatchParams
{
	enum BatchParamsVectors
	{
		INVALID = -1,
		ELITISM_PERC = 0,
		GENS_PER_WISDOM = 1,
		MAX_GENERATIONS = 2,
		MUTATION_PROBS = 3,
		CROSSOVER_PROBS = 4,
		POPULATION_SIZES = 5,
		MAX
	};

	string puzzleFilePath = "";
	vector<int> populationSizes;
	vector<float> crossoverProbs;
	vector<float> mutationProbs;
	vector<int> maxGenerations;
	vector<int> gensPerWisdoms;
	vector<float> elitismPercs;

	const int GetVectorSize(int index)
	{
		size_t size = 0;
		switch (static_cast<BatchParamsVectors>(index))
		{
		case BatchParamsVectors::ELITISM_PERC:
			size = elitismPercs.size();
			break;
		case BatchParamsVectors::GENS_PER_WISDOM:
			size = gensPerWisdoms.size();
			break;
		case BatchParamsVectors::MAX_GENERATIONS:
			size = maxGenerations.size();
			break;
		case BatchParamsVectors::MUTATION_PROBS:
			size = mutationProbs.size();
			break;
		case BatchParamsVectors::CROSSOVER_PROBS:
			size = crossoverProbs.size();
			break;
		case BatchParamsVectors::POPULATION_SIZES:
			size = populationSizes.size();
			break;
		default:
			return 0;
		}
		return static_cast<int>(size);
	}
};

struct BatchUpdateSave
{
	int currIndex = 0;
	vector<int> indices = { 0, 0, 0, 0, 0, 0 };
	int executionCount = 0;

	struct Parameters params;

	ofstream outputFile;
	bool bOutputHeaderUpdated = false;
};

/// <summary>
/// A struct used for Wisdom of Crowds.
/// </summary>
struct WoCConnection
{
	int islandID;
	bool bSingleConnection;
	int count;
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
	~Node();

	void PrintNodeInfo() const;

	int GetCurrentNumConnections() const;

	//A 'neighbor' is a struct that contains the # of bridges between the two nodes, the direction to the neighbor, as well as a pointer to the neighboring node object itself. 
	//Two nodes who are neighbors should have pointers to one another. 
	vector<Neighbor*> neighbors;

	//vector<Neighbor*> updatedNeighbors;
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
	//A vector of vectors that indicates the state of the board. Each element of the board dictates a different part of the board based on their value.
	//0:	An empty spacee on the board. 
	//1-8:	An island on the board, with a value as dictated with the number.
	//-1:	|	1 bridge vertical
	//-2:	||	2 bridges vertical
	//-3:	-	1 bridge horizontal
	//-4:	=	2 bridges horizontal
	vector<vector<int>> board;

	// Size of the board.
	int boardSizeX;
	int boardSizeY;

	vector<Node*> islands; //All the islands that are on the board.

public:
	/// <summary>
	/// Default constructor to help with initialization.
	/// </summary>
	HashiBoard() : bLongerWidth(false), boardSizeX(0), boardSizeY(0), currGen(0), bestPerc(0.f)
	{
		texture = SDL_CreateTexture(SDL->GetRenderer(), SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, SCREEN_WIDTH, SCREEN_HEIGHT);
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
	/// Initializes the hashi board with a file.
	/// </summary>
	/// <param name="filePath">Path to the board file.</param>
	/// <returns>Whether the board was initialized properly.</returns>
	bool Initialize(string filePath);

	/// <summary>
	/// Resets the current hashi board object.
	/// </summary>
	/// <returns>Whether the object was reset or not.</returns>
	bool Reset();

	/// <summary>
	/// Renders the hashi texture to the window.
	/// </summary>
	void RenderBoard();

	/// <summary>
	/// Updates the hashiboard by executing the algorithm.
	/// </summary>
	/// <param name="params">Parameters to pass to the algorithm.</param>
	/// <returns>Whether or not to continue updating.</returns>
	bool Update(Parameters params);

	/// <summary>
	/// Updates based off of batch testing.
	/// </summary>
	/// <param name="params">Batch parameters for batch testing.</param>
	/// <returns>Whether or not to continue executing.</returns>
	bool BatchUpdate(BatchParams params);

	/// <summary>
	/// Parses a string to a vector of integers.
	/// </summary>
	/// <param name="intsString">String to parse.</param>
	/// <param name="ints">Integers vector.</param>
	static void ParseIntString(string intsString, vector<int>& ints);

	/// <summary>
	/// Parses a string to a vector of floats.
	/// </summary>
	/// <param name="floatsString">String to parse.</param>
	/// <param name="floats">Floats vector.</param>
	static void ParseFloatString(string floatsString, vector<float>& floats);

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
	/// Current generation of the algorithm.
	/// </summary>
	int currGen;

	/// <summary>
	/// Best fit percentage of a chromosome.
	/// </summary>
	float bestPerc;

	/// <summary>
	/// Outpus csv of the algorithm.
	/// </summary>
	ofstream outputFile;

	/// <summary>
	/// Save data for batch testing.
	/// </summary>
	BatchUpdateSave batchSave;

	/// <summary>
	/// Start time for a single iteration.
	/// </summary>
	chrono::time_point<chrono::system_clock> startTime;

	/// <summary>
	/// Time of the current generation.
	/// </summary>
	chrono::time_point<chrono::system_clock> iterTime;

	// Useful typedef for less typing.
	typedef uint8_t uint8;
	// A gene is a pair of integer and uint8.
	// Integer defines the id of the island.
	// Uint8 defines the bridge connections:
	// Upper word - Defines double connections.
	// Lower word - Defines single connections.
	// E.x.
	// UDRL UDRL
	// 0100 1011
	// Upper word - Double connection down.
	// Lower word - Single connection up, right, and left.
	typedef pair<int, uint8> Gene;
	// A chromosome is a collection of genes and how close they are to a solution.
	typedef vector<Gene> Chromosome;
	typedef pair<float, Chromosome> FitnessChromosome;
	// A population is a collection of chromosomes.
	typedef vector<FitnessChromosome> Population;

	/// <summary>
	/// Population to be used in the algorithm.
	/// </summary>
	Population population;

private:
	/// <summary>
	/// Parses through a board file.
	/// </summary>
	/// <param name="filePath">Path to the board.</param>
	/// <returns>Whether or not the file was properly parsed.</returns>
	bool ParseBoardFile(string filePath);

	/// <summary>
	/// Parses through a given board in the board variable, and converts it to islands, and figures out the neighbors to the islands. 
	/// </summary>
	/// <returns>Returns true if board was successfully setup, false otherwise. </returns>
	bool Parseboard();

	/// <summary>
	/// Update a given node's neighbors by checking for the first 'hit' of a node in each cardinal direction.
	/// If a 'hit' gets nothing, or hits a bridge, then no neighbor exists in that direction. 
	/// </summary>
	/// <param name="node">Pointer to the node that is checking its neighbors.</param>
	/// <param name="bShouldClearNeighbors">Should the neighbors list be cleared out before updating? Should usually be true, but here just in case. </param>
	void UpdateBaseNeighborInfo(Node* node, bool bShouldClearNeighbors = true);

	void UpdateNeighborInfo(Node* node, Gene& gene, bool bShouldClearNeighbors = true);

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

	/// <summary>
	/// Processes the algorithm.
	/// </summary>
	/// <param name="params">Parameters to pass to the algorithm.</param>
	/// <returns>Whether or not to continue the algorithm.</returns>
	bool Process(Parameters params);

	/// <summary>
	/// Initializes the population for the algorithm.
	/// </summary>
	/// <param name="populationSize">Size of the requested population.</param>
	/// <param name="seed">Seed to see the random generator with.</param>
	/// <returns>Whether or not the population was created successfully.</returns>
	bool InitializePopulation(int populationSize, unsigned int seed);

	/// <summary>
	/// Initializes an island specified by the id parameter.
	/// </summary>
	/// <param name="id">Id of the island to initialize.</param>
	/// <param name="connection">It's bitmask to parse.</param>
	/// <param name="chrome">The chromosome it belongs to to duplicate connections to neighbors.</param>
	void InitializeIslandConnections(int id, uint8& connection, Chromosome& chrome);

	/// <summary>
	/// Repopulate the connections
	/// </summary>
	/// <param name="chrome">The chromosome to remake for</param>
	void RemakeIslandConnections(Chromosome& chrome);

	/// <summary>
	/// Fitness function.
	/// </summary>
	/// <param name="chrome">Chromosome to evaluate.</param>
	void EvaluateChromosome(FitnessChromosome& chrome);

	/// <summary>
	/// Calculates the number of connections from the bitmask.
	/// </summary>
	/// <param name="connection">Bitmask to parse.</param>
	/// <returns>The number of connections found in the bitmask.</returns>
	int CalcConnectionsFromMask(uint8 connection);

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

	/// <summary>
	/// Perform the crossover using single point crossover
	/// </summary>
	/// <param name="crossoverChromes"> The pair of </param>
	void PerformCrossover(const vector<pair<int, int>>& crossoverChromes);

	/// <summary>
	/// Perform the mutation: Flipping a bit randomly in the bridge connections
	/// Essentially adds or removes connections between two islands
	/// </summary>
	/// <param name="mutationChromes"> Chromosomes to mutate </param>
	void PerformMutation(const vector<int>& mutationChromes);

	/// <summary>
	/// Performs the Wisdom of Crowds operation on the population.
	/// </summary>
	/// <param name="elitismPerc">Percentage of the best of the population to consider.</param>
	void PerformWisdomOfCrowds(float elitismPerc);

	/// <summary>
	/// Clears the bridges on the board to prepare for rendering.
	/// </summary>
	void ClearBridgesOnBoard();

	/// <summary>
	/// Checks if the passed chromosome is unique to the population.
	/// </summary>
	/// <param name="chromosome">Chromosome to check.</param>
	/// <returns>Whether or not the chromosome is unique.</returns>
	bool CheckIfUnique(const Chromosome& chromosome);

	/// <summary>
	/// Fixes the chromosome connection errors that are present.
	/// </summary>
	/// <param name="chromosome">Chromosome to fix.</param>
	void FixChromosomeConnections(Chromosome& chromosome);

	/// <summary>
	/// Fixes the mirroring of direction connections to islands that are present.
	/// </summary>
	/// <param name="chromosome">Chromosome to fix.</param>
	void FixMirroringConnections(Chromosome& chromosome);

	void FixOverlappingConnections(Chromosome& chromosome);

	/// <summary>
	/// Fixes the excess connections that are present.
	/// </summary>
	/// <param name="chromosome">Chromosome to fix.</param>
	void FixExcessConnections(Chromosome& chromosome);

	//------------------------
	// none of these post processing methods seem to work ;-;

	void ReduceExcessConnections(uint8& mask, int excessConnections);

	void ValidateConnections(Chromosome& chrome);

	//------------------------

	/// <summary>
	/// Wisdom of crowds is complex, but essentially, it takes all of the 
	/// chromosomes at the end of a generations, and matches them up.
	/// </summary>
	void WisdomOfCrowds();

	/// <summary>
	/// Temporary Helper Function to print the board after every new generation 
	/// For debugging purposes - to be removed
	/// </summary>
	void PrintBoard() const;

	/// <summary>
	/// Return the opposite direction
	/// </summary>
	/// <param name="currentDirection"></param>
	/// <returns> Opposing direction </returns>
	Direction GetOppositeDirection(Direction currentDirection);

	/// <summary>
	/// See if the path is walkable
	/// </summary>
	/// <returns> True if path is separated, false otherwise </returns>
	bool isDisjoint() const;

	/// <summary>
	/// Converts a string to an integer.
	/// </summary>
	/// <param name="intString">String to convert.</param>
	/// <returns>Converted integer.</returns>
	static int StringToInt(string intString);

	/// <summary>
	/// Converts a string to a float.
	/// </summary>
	/// <param name="floatString">String to convert.</param>
	/// <returns>Converted float.</returns>
	static float StringToFloat(string floatString);
};

class WisdomOfCrowds
{
	HashiBoard board;
};