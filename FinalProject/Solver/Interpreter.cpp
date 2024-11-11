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

int Node::GetCurrentNumConnections() const
{
	int totalConnections = 0;
	for (Neighbor* neighbor : neighbors)
	{
		totalConnections += neighbor->numOfBridges;
	}
	return totalConnections;
}

Node::~Node()
{
	for (Neighbor* neighbor : neighbors)
	{
		delete neighbor;
	}
}

bool HashiBoard::Initialize(string filePath)
{
	if (ParseBoardFile(filePath))
	{
		bLongerWidth = board[0].size() >= board.size();

		boardSizeX = static_cast<int>(board[0].size());
		boardSizeY = static_cast<int>(board.size());

		Parseboard();
	}

	return false;
}

bool HashiBoard::Reset()
{
	if (!board.size())
	{
		return false;
	}

	SDL_DestroyTexture(texture);
	texture = SDL_CreateTexture(SDL->GetRenderer(), SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, SCREEN_WIDTH, SCREEN_HEIGHT);
	bLongerWidth = false;
	boardSizeX = 0;
	boardSizeY = 0;
	board.clear();
	islands.clear();
	population.clear();
	currGen = 0;
	bestPerc = 0.f;

	return true;
}

bool HashiBoard::ParseBoardFile(string filePath)
{
	// Open the file.
	ifstream file(filePath, ios::in);
	if (file.is_open())
	{
		// Helpful parsing variables.
		int numRows = 0;
		bool bBoardStarted = false;
		bool bRowStarted = false;

		// While the file hasn't ended.
		string line;
		while (!file.eof())
		{
			// Obtain a line and ensure it's valid.
			getline(file, line);
			if (line.size())
			{
				// Iterate through the line.
				string::const_iterator sIter = line.cbegin();
				for (sIter; sIter != line.cend(); ++sIter)
				{
					// If valid integer and in the board and row.
					int id = ASCII_ATOI(*sIter);
					if (id >= 0 && id <= 9 && bBoardStarted && bRowStarted)
					{
						board[numRows].push_back(id);
					}
					// If we see an open bracket and the board isn't started.
					else if (*sIter == '[' && !bBoardStarted)
					{
						bBoardStarted = true;
					}
					// If we see an open bracket and the board is started but the row isn't.
					else if (*sIter == '[' && bBoardStarted && !bRowStarted)
					{
						bRowStarted = true;
						board.push_back(vector<int>());
					}
					// If we see a closing bracket and the board is started and the row is started.
					else if (*sIter == ']' && bBoardStarted && bRowStarted)
					{
						++numRows;
						bRowStarted = false;
					}
					// If we see a closing bracket and the board is started but the row isn't.
					else if (*sIter == ']' && bBoardStarted && !bRowStarted)
					{
						bBoardStarted = false;
					}
				}
			}
		}

		file.close();
	}

	// If the board is populated, the parsing worked.
	if (board.size())
	{
		return true;
	}

	return false;
}

bool HashiBoard::Parseboard()
{
	int rowIndex = 0;
	int columnIndex = 0;
	int currentID = 0;
	for (const vector<int> row : board)
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
		UpdateBaseNeighborInfo(checkedNode);
		checkedNode->PrintNodeInfo();
	}
	return true;
}

void HashiBoard::UpdateBaseNeighborInfo(Node* node, bool bShouldClearNeighbors)
{;
	Node* neighborNode = nullptr;
	if (bShouldClearNeighbors)
	{
		for (Neighbor* neighbor : node->neighbors) 
		{
			delete neighbor; 
		}
		node->neighbors.clear();
	}
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
		for (index = col - 1; index >= 0; index--)
		{
			parsedNodeValue = board[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue != HORI_SINGLE_BRIDGE && parsedNodeValue != HORI_DOUBLE_BRIDGE && parsedNodeValue != 0) { return nullptr; } //A vertical bridge has been hit. 
		}
		break;
	case Direction::RIGHT:
		for (index = col + 1; index < boardSizeX; index++)
		{
			parsedNodeValue = board[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue != HORI_SINGLE_BRIDGE && parsedNodeValue != HORI_DOUBLE_BRIDGE && parsedNodeValue != 0) { return nullptr; } //A vertical bridge has been hit. 
		}
		break;
	case Direction::DOWN:
		for (index = row + 1; index < boardSizeY; index++)
		{
			parsedNodeValue = board[index][col];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue != VERT_SINGLE_BRIDGE && parsedNodeValue != VERT_DOUBLE_BRIDGE && parsedNodeValue != 0) { return nullptr; } //A horizontal bridge has been hit. 
		}
		break;
	case Direction::UP:
		for (index = row - 1; index >= 0; index--)
		{
			parsedNodeValue = board[index][col];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue != VERT_SINGLE_BRIDGE && parsedNodeValue != VERT_DOUBLE_BRIDGE && parsedNodeValue != 0) { return nullptr; } //A horizontal bridge has been hit. 
		}
		break;
	}
	return nullptr;
}

void HashiBoard::ClearBridgesOnBoard()
{
	for (int y = 0; y < board.size(); ++y)
	{
		for (int x = 0; x < board[0].size(); ++x)
		{
			if (board[y][x] < 0)
			{
				board[y][x] = 0;
			}
		}
	}
}

bool HashiBoard::BatchUpdate(BatchParams params)
{
	// Initialization of the batch.
	if (!batchSave.params.seed)
	{
		batchSave.params = { static_cast<unsigned int>(time(0)), 
			params.populationSizes[0], 
			params.crossoverProbs[0], 
			params.mutationProbs[0], 
			params.maxGenerations[0],
			false, 
			params.gensPerWisdoms[0], 
			params.elitismPercs[0]};

		batchSave.outputFile.open("BatchOutput.csv");
		if (batchSave.outputFile.is_open())
		{
			batchSave.outputFile << "Execution Count,Seed,Generation,Solved Perc,Sec Start,Sec Iter," << 
				"Population Size,Crossover Prob,Mutation Prob,Max Generations,With Wisdom,Gens Per Wisdom,Elitism Perc" << endl;
		}
		else
		{
			cout << "Error Opening Batch Output File!" << endl;
			return false;
		}
	}
	else
	{
		batchSave.params = { static_cast<unsigned int>(time(0)),
			params.populationSizes[batchSave.indices[BatchParams::BatchParamsVectors::POPULATION_SIZES]],
			params.crossoverProbs[batchSave.indices[BatchParams::BatchParamsVectors::CROSSOVER_PROBS]],
			params.mutationProbs[batchSave.indices[BatchParams::BatchParamsVectors::MUTATION_PROBS]],
			params.maxGenerations[batchSave.indices[BatchParams::BatchParamsVectors::MAX_GENERATIONS]],
			(batchSave.executionCount % 2 != 0 ? true : false),
			params.gensPerWisdoms[batchSave.indices[BatchParams::BatchParamsVectors::GENS_PER_WISDOM]],
			params.elitismPercs[batchSave.indices[BatchParams::BatchParamsVectors::ELITISM_PERC]] };
	}

	// If we need to reprint the header data to the output file.
	if (!batchSave.bOutputHeaderUpdated)
	{
		Parameters& refParams = batchSave.params;
		batchSave.outputFile << ++batchSave.executionCount << ',' << refParams.seed << ",,,,," <<
			refParams.populationSize << ',' << refParams.crossoverProb << ',' << refParams.mutationProb << ',' << refParams.maxGenerations << ',' <<
			refParams.bWithWisdom << ',' << refParams.gensPerWisdom << ',' << refParams.elitismPerc << endl;
		batchSave.bOutputHeaderUpdated = true;

		// Helpful visualization in the console.
		cout << "----- BATCH TESTING -----" << endl
			<< "Execution Count:\t" << batchSave.executionCount << endl
			<< "Population Size:\t" << batchSave.params.populationSize << endl
			<< "Crossover Probability:\t" << batchSave.params.crossoverProb << endl
			<< "Mutation Probability:\t" << batchSave.params.mutationProb << endl
			<< "Maximum Generations:\t" << batchSave.params.maxGenerations << endl
			<< "With Wisdom:\t\t" << batchSave.params.bWithWisdom << endl
			<< "Generations per Wisdom:\t" << batchSave.params.gensPerWisdom << endl
			<< "Elitism Percentage:\t" << batchSave.params.elitismPerc << endl
			<< "-------------------------" << endl;
	}

	// Base update to just run 1 execution at a time.
	bool bContinue = Update(batchSave.params);

	// I fthe execution has ended.
	if (!bContinue)
	{
		// Obtain the total time and output to the file.
		chrono::duration<double> elapsedSinceStart = chrono::system_clock::now() - startTime;
		batchSave.outputFile << "Total Time," << elapsedSinceStart.count() << endl;

		// Reset the header data.
		batchSave.bOutputHeaderUpdated = false;
		// If we haven't ran on wisdom yet, ignore.
		if (batchSave.executionCount % 2 == 0)
		{
			// Increment the current parameters testing.
			if (batchSave.indices[batchSave.currIndex] + 1 < params.GetVectorSize(batchSave.currIndex))
			{
				++batchSave.indices[batchSave.currIndex];
			}
			// If we can't, move to the next parameter.
			else
			{
				// Iterate till a parameter that can be changed is found.
				++batchSave.currIndex;
				while (batchSave.currIndex < static_cast<int>(BatchParams::BatchParamsVectors::MAX)
					&& batchSave.indices[batchSave.currIndex] + 1 >= params.GetVectorSize(batchSave.currIndex))
				{
					++batchSave.currIndex;
				}
				// If found.
				if (batchSave.currIndex < static_cast<int>(BatchParams::BatchParamsVectors::MAX))
				{
					// Reset all previous parameters and continue.
					++batchSave.indices[batchSave.currIndex];
					for (--batchSave.currIndex; batchSave.currIndex >= 0; --batchSave.currIndex)
					{
						batchSave.indices[batchSave.currIndex] = 0;
					}
					++batchSave.currIndex;
				}
			}
		}

		// If a parameter is found.
		if (batchSave.currIndex < static_cast<int>(BatchParams::BatchParamsVectors::MAX))
		{
			// Reset, initialize, and continue executing.
			Reset();
			Initialize(params.puzzleFilePath);
			bContinue = true;
		}
		else
		{
			// Otherwise, close the file and do nothing to bContinue.
			batchSave.outputFile.close();
			batchSave.bOutputHeaderUpdated = false;
			batchSave.currIndex = 0;
			batchSave.executionCount = 0;
			for (int& index : batchSave.indices)
			{
				index = 0;
			}
			batchSave.params = Parameters();
		}
	}

	return bContinue;
}

bool HashiBoard::Update(Parameters params)
{
	// Process the algorithm.
	bool bContinue = Process(params);

	// Find the chromosome with the best completion percentage.
	Population::iterator Iter = find_if(population.begin(), population.end(),
		[this](const FitnessChromosome& a)
		{
			return a.first == bestPerc;
		});
	Chromosome& chrome = (*Iter).second;

	ClearBridgesOnBoard();

	// Iterate through the genes to populate the board.
	for (Gene& gene : chrome)
	{
		uint8 mask = gene.second;
		Node* node = islands[gene.first];

		int x = node->coords[1];
		int y = node->coords[0];

		int bitCount = 0;
		while (mask)
		{
			if (mask & 1)
			{
				switch (bitCount % BITMASK_BOUNDARY)
				{
				case 0: // Left
					for (int index = x - 1; board[y][index] <= 0; --index)
					{
						board[y][index] = (bitCount < BITMASK_BOUNDARY) ? HORI_SINGLE_BRIDGE : HORI_DOUBLE_BRIDGE;
					}
					break;
				case 1: // Right
					for (int index = x + 1; board[y][index] <= 0; ++index)
					{
						board[y][index] = (bitCount < BITMASK_BOUNDARY) ? HORI_SINGLE_BRIDGE : HORI_DOUBLE_BRIDGE;
					}
					break;
				case 2: // Down
					for (int index = y + 1; board[index][x] <= 0; ++index)
					{
						board[index][x] = (bitCount < BITMASK_BOUNDARY) ? VERT_SINGLE_BRIDGE : VERT_DOUBLE_BRIDGE;
					}
					break;
				case 3: // Up
					for (int index = y - 1; board[index][x] <= 0; --index)
					{
						board[index][x] = (bitCount < BITMASK_BOUNDARY) ? VERT_SINGLE_BRIDGE : VERT_DOUBLE_BRIDGE;
					}
					break;
				default:
					break;
				}
			}
			mask >>= 1;
			++bitCount;
		}
	}

	return bContinue;
}

bool HashiBoard::Process(Parameters params)
{
	// Initialization.
	if (!population.size())
	{
		InitializePopulation(params.populationSize, params.seed);
	}
	++currGen;
	vector<pair<int, int>> crossoverChromes;
	bool bCrossoverStarted = false;
	vector<int> mutationChromes;
	// Determining crossovers and mutations.
	for (int i = 0; i < population.size(); ++i)
	{
		float perc = static_cast<float>(rand()) / RAND_MAX;
		if (perc < params.crossoverProb/* && population[i].first >= 0.6f*/)
		{
			if (!bCrossoverStarted)
			{
				bCrossoverStarted = true;
				crossoverChromes.push_back(pair<int, int>(i, -1));
			}
			else
			{
				bCrossoverStarted = false;
				crossoverChromes[crossoverChromes.size() - 1].second = i;
			}
		}
		if (perc < params.mutationProb)
		{
			mutationChromes.push_back(i);
		}
	}
	// Performing crossovers and mutations.
	if (crossoverChromes.size())
	{
		PerformCrossover(crossoverChromes);
	}

	if (mutationChromes.size())
	{
		PerformMutation(mutationChromes);
	}
	// Evaluate fitness.
	for (FitnessChromosome& chrome : population)
	{
		EvaluateChromosome(chrome);
	}
	// Removing bad parents.
	sort(population.begin(), population.end(),
		[this](const FitnessChromosome& a, const FitnessChromosome& b)
		{
			return a.first > b.first;
		});
	if (population.size() > params.populationSize)
	{
		population.erase(population.begin() + params.populationSize, population.end());
	}
	// Perform Wisdom of Crowds if enabled.
	if (params.bWithWisdom && currGen % params.gensPerWisdom == 0)
	{
		PerformWisdomOfCrowds(params.elitismPerc);

		sort(population.begin(), population.end(),
			[this](const FitnessChromosome& a, const FitnessChromosome& b)
			{
				return a.first > b.first;
			});
	}
	// Outputting to csv if new best parent.
	if (population[0].first != bestPerc)
	{
		bestPerc = population[0].first;
		// Obtain some time data.
		chrono::time_point<chrono::system_clock> now = chrono::system_clock::now();
		chrono::duration<double> elapsedSinceStart = now - startTime;
		chrono::duration<double> elapsedSinceIter = now - iterTime;
		iterTime = now;
		// If not batch testing.
		if(!batchSave.outputFile.is_open())
		{
			outputFile << currGen << ',' << bestPerc << ',' << elapsedSinceStart.count() << ',' << elapsedSinceIter.count() << endl;
		}
		else
		{
			batchSave.outputFile << ",," << currGen << ',' << bestPerc << ',' << elapsedSinceStart.count() << ',' << elapsedSinceIter.count() << ",,,,,,, " << endl;
		}
	}
	// Ending condition.
	if (bestPerc >= 1.f || currGen >= params.maxGenerations)
	{
		chrono::duration<double> elapsedSinceStart = chrono::system_clock::now() - startTime;
		// If not batch testing.
		if(!batchSave.outputFile.is_open())
		{
			outputFile << "Total Time," << elapsedSinceStart.count() << endl;
			outputFile.close();
		}
		
		cout << "Generation: " << currGen << "| Fitness: " << bestPerc << endl;
		PrintBoard();
		return false;
	}

	return true;
}

bool HashiBoard::InitializePopulation(int populationSize, unsigned int seed)
{
	// Seed the random generator.
	cout << "Seed: " << seed << endl;
	srand(seed);
	// Obtain time data.
	startTime = chrono::system_clock::now();
	iterTime = startTime;
	// If not batch testing.
	if(!batchSave.outputFile.is_open())
	{
		outputFile.open("Output.csv");
		if (outputFile.is_open())
		{
			outputFile << "Seed," << seed << endl;
			outputFile << "Generation,Solved Perc,Sec Start,Sec Iter" << endl;
		}
	}

	// Create an empty chromosome.
	Chromosome emptyChrome;
	for (Node* node : islands)
	{
		emptyChrome.push_back(Gene(node->nodeID, 0));
	}

	// Iterate till the population is filled.
	int count = 0;
	while (count < populationSize)
	{
		// Add an empty chromosome.
		population.push_back(FitnessChromosome(0.f, Chromosome(emptyChrome)));
		Chromosome& chromosome = population[count].second;
		//Iterate through the genes and initialize connections.
		for (int i = 0; i < chromosome.size(); ++i)
		{
			Gene& gene = chromosome[i];
			InitializeIslandConnections(gene.first, gene.second, chromosome);
		}
		FixChromosomeConnections(population[count].second);
		//Evaluate the fitness of the chromosome before leaving.
		EvaluateChromosome(population[count]);
		++count;
	}

	return true;
}

void HashiBoard::InitializeIslandConnections(int id, uint8& connection, Chromosome& chrome)
{
	// Obtain the island and some helpful values.
	Node* island = islands[id];
	int numConnections = CalcConnectionsFromMask(connection);
	int numIterations = static_cast<int>(islands.size());
	// Iterate till the max number of iterations or the island has reached max connections.
	while (numIterations-- && numConnections < island->value)
	{
		// Obtain a random neighbor and associated values.
		if(island->neighbors.size() <= 0) continue;
		int randNeighbor = static_cast<int>(rand() % island->neighbors.size());
		Neighbor* nb = island->neighbors[randNeighbor];
		Chromosome::iterator Iter = find_if(chrome.begin(), chrome.end(),
			[nb](const Gene& a)
			{
				return a.first == nb->neighborNode->nodeID;
			});
		Gene* nbGene = nullptr;
		if (Iter != chrome.end())
		{
			nbGene = &(*Iter);
		}
		if (nbGene)
		{
			// Obtain the open connections.
			int remainingConnections = island->value - numConnections;
			if (remainingConnections > 0)
			{
				// Calculate how many connections from the current island to the neighbor there can be.
				int valDiff = nb->neighborNode->value - CalcConnectionsFromMask(nbGene->second) - numConnections;
				if (valDiff >= 2)
				{
					if (remainingConnections >= 2)
					{
						// If there's already a double bridge in this direction.
						if (connection & (1 << static_cast<int>(nb->neighborDirection) << BITMASK_BOUNDARY))
						{
							continue;
						}
						// If not, add one.
						connection |= (1 << static_cast<int>(nb->neighborDirection) << BITMASK_BOUNDARY);
						//nb->numOfBridges = 2;
						numConnections += 2;
						nbGene->second |= (1 << (static_cast<int>(nb->neighborDirection) ^ 1) << BITMASK_BOUNDARY);
					}
					else
					{
						// If there's already a single bridge in this direction.
						if (connection & 1 << static_cast<int>(nb->neighborDirection))
						{
							continue;
						}
						// If not, add one.
						connection |= 1 << static_cast<int>(nb->neighborDirection);
						//nb->numOfBridges = 1;
						++numConnections;
						nbGene->second |= 1 << (static_cast<int>(nb->neighborDirection) ^ 1);
					}
				}
				else if (valDiff > 0)
				{
					// If there's already a single bridge in this direction.
					if (connection & 1 << static_cast<int>(nb->neighborDirection))
					{
						continue;
					}
					// If not, add one.
					connection |= 1 << static_cast<int>(nb->neighborDirection);
					//nb->numOfBridges = 1;
					++numConnections;
					nbGene->second |= 1 << (static_cast<int>(nb->neighborDirection) ^ 1);
				}
			}
		}
	}
}

void HashiBoard::RemakeIslandConnections(Chromosome& chrome)
{
	// Obtain the island and some helpful values.
	ClearBridgesOnBoard();

	Chromosome shuffledChrome = chrome;
	random_shuffle(shuffledChrome.begin(), shuffledChrome.end());

	// Iterate through the genes to populate the board.
	for (Gene& gene : shuffledChrome)
	{
		uint8 mask = gene.second;
		Node* node = islands[gene.first];

		int x = node->coords[1];
		int y = node->coords[0];

		int bitCount = 0;
		while (mask)
		{
			if (mask & 1)
			{
				switch (bitCount % BITMASK_BOUNDARY)
				{
				case 0: // Left
					for (int index = x - 1; board[y][index] <= 0; --index)
					{
						board[y][index] = (bitCount < BITMASK_BOUNDARY) ? HORI_SINGLE_BRIDGE : HORI_DOUBLE_BRIDGE;
					}
					break;
				case 1: // Right
					for (int index = x + 1; board[y][index] <= 0; ++index)
					{
						board[y][index] = (bitCount < BITMASK_BOUNDARY) ? HORI_SINGLE_BRIDGE : HORI_DOUBLE_BRIDGE;
					}
					break;
				case 2: // Down
					for (int index = y + 1; board[index][x] <= 0; ++index)
					{
						board[index][x] = (bitCount < BITMASK_BOUNDARY) ? VERT_SINGLE_BRIDGE : VERT_DOUBLE_BRIDGE;
					}
					break;
				case 3: // Up
					for (int index = y - 1; board[index][x] <= 0; --index)
					{
						board[index][x] = (bitCount < BITMASK_BOUNDARY) ? VERT_SINGLE_BRIDGE : VERT_DOUBLE_BRIDGE;
					}
					break;
				default:
					break;
				}
			}
			mask >>= 1;
			++bitCount;
		}
	}
}

void HashiBoard::EvaluateChromosome(FitnessChromosome& chrome)
{
	// Fitness is defined as the completion of all island connections, as well as penalizing 'wrong' connections. 
	int numRequiredConnections = 0;
	int numAcquiredConnections = 0;

	//Remake the board using the specifications of the chromesome. 
	RemakeIslandConnections(chrome.second);

	//We give a harsh penalty to any nodes that are in an 'impossible' state.
	int impossibleNodes = 0;
	float penaltyImpossibleNode = -0.05f;

	//We have a penalty for any bridges that were overlapping.
	int overlappingBridges = 0;
	float penaltyOverlappingBridge = -0.02f;

	//We have a small penalty for any nodes that have not reached their target value. 
	int incorrectNodeValues = 0;
	float penaltyExcessiveNode = -0.01f;
	float penaltyDisjoint = -0.1f;
	//Now, we need to update each island
	for (Node* island : islands)
	{
		Gene& gene = chrome.second[island->nodeID];
		Node* neighborNode = nullptr;
		for (Neighbor* neighbor : island->neighbors)
		{
			delete neighbor;
		}
		island->neighbors.clear();

		//Now we check if any connections were severed compared to the gene. If so, then this is due to an overlap! 
		uint8 mask = gene.second;
		//We have an array of bools to signify if there is a connection in a specific direction.
		//The index for each direction in geneConnections is below:
		//Left = 0
		//Right = 1
		//Down = 2
		//Up = 3;
		int numGeneConnections[4] = { 0,0,0,0 };
		bool geneConnections[4] = { 0,0,0,0 };
		int bitCount = 0;
		while (mask)
		{
			if (mask & 1)
			{
				geneConnections[bitCount % BITMASK_BOUNDARY] = true;
				numGeneConnections[bitCount % BITMASK_BOUNDARY] = (bitCount < BITMASK_BOUNDARY) ? 1 : 2;
			}
			mask >>= 1;
			++bitCount;
		}

		//We need to traverse in each direction in the board in order to find a hit.
		Direction directionsToCheck[4] = { Direction::LEFT, Direction::RIGHT, Direction::DOWN, Direction::UP };
		for (Direction direction : directionsToCheck)
		{
			neighborNode = GetNodeInDirection(direction, island->coords[0], island->coords[1]);
			if (neighborNode)
			{
				island->neighbors.push_back(new Neighbor(0, direction, neighborNode));
				island->neighbors.back()->numOfBridges = numGeneConnections[(int)direction];

				//Since we care about properly making sure connections are made, make the neighbor have this node as its neighbor too, as long as they're not already neighbors.
				bool bAlreadyNeighbors = false;
				//Check the neighbor's neighbors, see if this node is already included. 
				for (Neighbor* neighborsNeighbor : neighborNode->neighbors)
				{
					if (neighborsNeighbor->neighborNode->nodeID == island->nodeID)
					{
						//If so, just validate the connection, and make sure the number of bridges between them is updated! 
						bAlreadyNeighbors = true;
						if (neighborsNeighbor->numOfBridges != island->neighbors.back()->numOfBridges)
						{
							neighborsNeighbor->numOfBridges = island->neighbors.back()->numOfBridges;
						}
					}
				}
				if (!bAlreadyNeighbors)
				{
					neighborNode->neighbors.push_back(new Neighbor(0, GetOppositeDirection(direction), island));
					neighborNode->neighbors.back()->numOfBridges = numGeneConnections[(int)direction];
				}

			}
		}

		//Check if an island has too many or too few nodes on it; give it a penalty if it does.
		if (island->GetCurrentNumConnections() != island->value)
		{
			incorrectNodeValues++;
		}

		//Check for cases of nodes with a value of 1 or 2 attaching to another node of equal value.
		//That's an impossible case to have! Penalize that as well.
		if (island->value <= 2)
		{
			for (Neighbor* neighbor : island->neighbors)
			{
				if (neighbor->numOfBridges > 0)
				{
					if (island->value == neighbor->neighborNode->value && island->value == neighbor->numOfBridges)
					{
						impossibleNodes++;
					}
				}
			}
		}

		//Now, we iterate through each neighbor in the island; for a neighbor that corresponds to
		//a direction in the geneConnections array, we set that direction in the geneConnections array to be false.
		for (Neighbor* neighbor : island->neighbors)
		{
			geneConnections[(int)neighbor->neighborDirection] = 0;
		}

		//After this process is done, if there are still any 'true' values in the geneConnections array,
		//then we know we've lost a connection in there due to an overlapping bridge.
		for (bool lostConnection : geneConnections)
		{
			if (lostConnection == true)
			{
				overlappingBridges++;
			}
		}

		//Check if an island is incomplete; can it stil be completed in this new board state?
		if (island->value > island->GetCurrentNumConnections())
		{
			//The island has no more neighbors! That's a no-no. 
			if (island->neighbors.size() == 0)
			{
				impossibleNodes++;
			}
			//The island doesn't have enough neighbors to fufill its value!
			//For example, a '3' island with only 1 neighbor could only fill 2 bridges. 
			if (island->neighbors.size() * 2 < island->value)
			{
				impossibleNodes++;
			}
		}
	}

	for (Gene& gene : chrome.second)
	{
		numRequiredConnections += islands[gene.first]->value;
		numAcquiredConnections += CalcConnectionsFromMask(gene.second);
	}

	int disjointCount = IsDisjoint();

	chrome.first = static_cast<float>(numAcquiredConnections) / numRequiredConnections;
	chrome.first += penaltyImpossibleNode * impossibleNodes;
	chrome.first += penaltyOverlappingBridge * overlappingBridges;
	chrome.first += penaltyExcessiveNode * incorrectNodeValues;
	chrome.first += penaltyDisjoint * disjointCount;
}

int HashiBoard::CalcConnectionsFromMask(uint8 connection)
{
	// Obtain the mask and other values.
	uint8 mask = connection;
	int bitCount = 0;
	int count = 0;
	// While the mask is nonzero.
	while (mask)
	{
		// If we are in double bridge territory.
		if (bitCount < BITMASK_BOUNDARY)
		{
			count += mask & 1;
		}
		// If we are in single bridge territory.
		else
		{
			count += (mask & 1) * 2;
		}
		mask >>= 1;
		++bitCount;
	}

	return count;
}

void HashiBoard::RenderBoard()
{
	if (!board.size())
	{
		return;
	}

	// Swap render target to the member texture.
	SDL_SetRenderTarget(SDL->GetRenderer(), texture);
	// Clear the texture.
	SDL_RenderClear(SDL->GetRenderer());
	// Set the color to white.
	SDL_SetRenderDrawColor(SDL->GetRenderer(), 255, 255, 255, 255);

	// Obtain the grid width and height.
	const int gridHeight = GRID_HEIGHT(board.size());
	const int gridWidth = GRID_WIDTH(board[0].size());

	// Iterate over the game board.
	for (int y = 0; y < board.size(); ++y)
	{
		for (int x = 0; x < board[0].size(); ++x)
		{
			int gridId = board[y][x];

			// If the id is positive, an island is formed.
			if (gridId > 0)
			{
				RenderIsland(gridId, x * gridWidth + gridWidth / 2, y * gridHeight + gridHeight / 2);
			}
			// If the id is negative, a bridge is formed.
			else if (gridId < 0)
			{
				RenderBridge(gridId, x * gridWidth + gridWidth / 2, y * gridHeight + gridHeight / 2);
			}
		}
	}

	// Change back to the default render target.
	SDL_SetRenderTarget(SDL->GetRenderer(), NULL);
	// Set the color back to black.
	SDL_SetRenderDrawColor(SDL->GetRenderer(), 0, 0, 0, 255);
	// Copy the texture.
	SDL_RenderCopy(SDL->GetRenderer(), texture, NULL, NULL);
}

const int HashiBoard::GetIslandRadius() const
{
	return (bLongerWidth ? GRID_WIDTH(board[0].size()) : GRID_HEIGHT(board.size())) * ISLAND_RADIUS_FACTOR;
}

void HashiBoard::RenderIsland(const int Id, const int CenterX, const int CenterY)
{
	// Obtains the radius and radius squared.
	const int r = GetIslandRadius();
	const double r2 = pow(r, 2);

	// Populates a vector with SDL_Points that are on the radius of the island.
	vector<SDL_Point> points;
	for (int y = CenterY - r; y < CenterY + r + 1; ++y)
	{
		for (int x = CenterX - r; x < CenterX + r + 1; ++x)
		{
			double result = pow(x - CenterX, 2) + pow(y - CenterY, 2);
			if (r2 - r <= result && result <= r2 + r)
			{
				SDL_Point point = { x, y };
				points.push_back(point);
			}
		}
	}
	// Draws the points as the circle.
	SDL_RenderDrawPoints(SDL->GetRenderer(), points.data(), static_cast<int>(points.size()));

	// Creates a null-terminated string of the id. NOTE: Ids don't go above 2 digits so this works.
	vector<char> idString = { ASCII_ITOA(Id), '\0' };
	// Creates a surfaces from the id and font.
	SDL_Surface* TextSurface = TTF_RenderText_Solid(SDL->GetFont(), idString.data(), { 255, 255, 255, 255 });
	if (TextSurface)
	{
		// Creates a texture from the surface.
		SDL_Texture* TextTexture = SDL_CreateTextureFromSurface(SDL->GetRenderer(), TextSurface);
		if (TextTexture)
		{
			// Calculates the render quad.
			SDL_Rect renderQuad = { CenterX - r * OFFSET_FACTOR, CenterY - r * OFFSET_FACTOR, TextSurface->w, TextSurface->h };
			// Renders the texture to the board texture.
			SDL_RenderCopy(SDL->GetRenderer(), TextTexture, NULL, &renderQuad);

			// Destroys references for no memory bloating.
			SDL_FreeSurface(TextSurface);
			SDL_DestroyTexture(TextTexture);
		}
	}
}

void HashiBoard::RenderBridge(const int Type, const int CenterX, const int CenterY)
{
	// Obtains the radius of the islands.
	const int r = GetIslandRadius();
	// Obtains the grid width and height.
	const int gridHeight = GRID_HEIGHT(board.size());
	const int gridWidth = GRID_WIDTH(board[0].size());
	// Calculates the shift in both directions needed.
	const int horizontalIslandShift = bLongerWidth ? r * 2 : static_cast<int>(r * BRIDGE_LENGTH_FACTOR);
	const int verticalIslandShift = !bLongerWidth ? r * 2 : static_cast<int>(r * BRIDGE_LENGTH_FACTOR);

	// Different logic depending on the type of bridge.
	switch (Type)
	{
	case VERT_SINGLE_BRIDGE: // |
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX, CenterY - gridHeight / 2 - verticalIslandShift,
			CenterX, CenterY + gridHeight / 2 + verticalIslandShift);
		break;
	case VERT_DOUBLE_BRIDGE: // ||
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX - r * OFFSET_FACTOR, CenterY - gridHeight / 2 - verticalIslandShift,
			CenterX - r * OFFSET_FACTOR, CenterY + gridHeight / 2 + verticalIslandShift);
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX + r * OFFSET_FACTOR, CenterY - gridHeight / 2 - verticalIslandShift,
			CenterX + r * OFFSET_FACTOR, CenterY + gridHeight / 2 + verticalIslandShift);
		break;
	case HORI_SINGLE_BRIDGE: // -
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX - gridWidth / 2 - horizontalIslandShift, CenterY,
			CenterX + gridWidth / 2 + horizontalIslandShift, CenterY);
		break;
	case HORI_DOUBLE_BRIDGE: // =
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX - gridWidth / 2 - horizontalIslandShift, CenterY - r * OFFSET_FACTOR,
			CenterX + gridWidth / 2 + horizontalIslandShift, CenterY - r * OFFSET_FACTOR);
		SDL_RenderDrawLine(SDL->GetRenderer(),
			CenterX - gridWidth / 2 - horizontalIslandShift, CenterY + r * OFFSET_FACTOR,
			CenterX + gridWidth / 2 + horizontalIslandShift, CenterY + r * OFFSET_FACTOR);
		break;
	default:
		break;
	}
}

void HashiBoard::PerformCrossover(const vector<pair<int, int>>& crossoverChromes)
{
	for (const auto& pair : crossoverChromes)
	{
		if (pair.second == -1) continue;

		// Get parent chromosomes
		Chromosome& parent1 = population[pair.first].second;
		Chromosome& parent2 = population[pair.second].second;

		// Create two child chromosomes by crossover
		Chromosome child1, child2;
		int numGenes = static_cast<int>(parent1.size());
		int crossoverPoint = rand() % numGenes;

		// First child (parent1's genes before crossover point, then parent2's)
		for (int i = 0; i < numGenes; ++i) {
			if (i < crossoverPoint) {
				child1.push_back(parent1[i]);
				child2.push_back(parent2[i]);
			}
			else {
				child1.push_back(parent2[i]);
				child2.push_back(parent1[i]);
			}
		}

		// Fix chromosomes to ensure they respect constraints
		FixChromosomeConnections(child1);
		FixChromosomeConnections(child2);

		if (CheckIfUnique(child1))
		{
			FitnessChromosome fChild1 = { 0.f, child1 };
			EvaluateChromosome(fChild1);
			population.push_back(fChild1);
		}

		if (CheckIfUnique(child2))
		{
			FitnessChromosome fChild2 = { 0.f, child2 };
			EvaluateChromosome(fChild2);
			population.push_back(fChild2);
		}
	}
}

vector<int> HashiBoard::GetUnsatisfiedNodes(const Chromosome& chromosome) {
	vector<int> unsatisfiedNodes;
	for (int i = 0; i < islands.size(); ++i) {
		Node* node = islands[i];
		int totalConnections = CalcConnectionsFromMask(chromosome[i].second);
		if (totalConnections != node->value) {
			unsatisfiedNodes.push_back(i);
		}
	}
	return unsatisfiedNodes;
}
void HashiBoard::PerformMutation(const vector<int>& mutationChromes) {
	for (int chromeIndex : mutationChromes) {
		Chromosome& chromosome = population[chromeIndex].second;
		vector<int> unsatisfiedNodes = GetUnsatisfiedNodes(chromosome);

		// Clear bridges and update base neighbor info
		ClearBridgesOnBoard();
		for (Node* node : islands) {
			UpdateBaseNeighborInfo(node);
		}

		// Only mutate unsatisfied nodes to attempt valid connections
		for (int nodeId : unsatisfiedNodes) {
			Node* node = islands[nodeId];
			Gene& gene = *find_if(chromosome.begin(), chromosome.end(),
				[nodeId](const Gene& g) { return g.first == nodeId; });

			// Toggle connections for the unsatisfied node
			uint8& connection = gene.second;
			int bitToToggle;
			Neighbor* nb = nullptr;
			bool bFoundBit = false;
			do {
				bitToToggle = rand() % (2 * BITMASK_BOUNDARY);
				Direction dir = static_cast<Direction>(bitToToggle % BITMASK_BOUNDARY);

				auto iter = find_if(node->neighbors.begin(), node->neighbors.end(),
					[&dir](const Neighbor* n) { return n->neighborDirection == dir; });
				if (iter == node->neighbors.end()) {
					bFoundBit = false;
				}
				else {
					bFoundBit = true;
					nb = *iter;
				}
			} while (!bFoundBit && !nb);

			// Ensure we update the bridge correctly
			int singleBit = bitToToggle % BITMASK_BOUNDARY;
			int doubleBit = singleBit + BITMASK_BOUNDARY;
			if ((connection & (1 << singleBit)) && !(connection & (1 << doubleBit))) {
				// If there's a single bridge, upgrade to a double bridge
				connection |= (1 << doubleBit);  // Set double bridge bit
				connection &= ~(1 << singleBit); // Clear single bridge bit
			}
			else if (!(connection & (1 << singleBit)) && (connection & (1 << doubleBit))) {
				// If there's a double bridge, downgrade to a single bridge
				connection &= ~(1 << doubleBit); // Clear double bridge bit
				connection |= (1 << singleBit);  // Set single bridge bit
			}
			else {
				// Toggle the selected bit as per mutation logic
				connection ^= (1 << bitToToggle);
			}

			// Update the neighbor’s gene for consistency, ensuring alignment
			if (nb) {
				int neighborShift = (bitToToggle < BITMASK_BOUNDARY ?
					((bitToToggle % BITMASK_BOUNDARY) ^ 1) :
					(((bitToToggle % BITMASK_BOUNDARY) ^ 1) + BITMASK_BOUNDARY));
				Gene* nbGene = &(*find_if(chromosome.begin(), chromosome.end(),
					[&nb](const Gene& g) { return g.first == nb->neighborNode->nodeID; }));

				// Align neighbor bit to ensure symmetric connection
				if (connection & (1 << bitToToggle)) {
					nbGene->second |= (1 << neighborShift);  // Set the bit on the neighbor
				}
				else {
					nbGene->second &= ~(1 << neighborShift); // Clear the bit on the neighbor
				}
			}
		}

		// Fix connections to ensure all bridge toggles are consistent
		FixChromosomeConnections(chromosome);
		EvaluateChromosome(population[chromeIndex]);
	}
}



void HashiBoard::PerformWisdomOfCrowds(float elitismPerc)
{
	// Map to hold the amount of connections per node.
	map<int, vector<WoCConnection>> connectionsMap;
	
	// Iterate through each chromosome.
	int count = 0;
	for (FitnessChromosome& fchrome : population)
	{
		if (count > population.size() * elitismPerc)
		{
			break;
		}

		Chromosome& chrome = fchrome.second;
		for (int i = 1; i < chrome.size(); ++i)
		{
			// Obtain the gene and the connection.
			const Gene& gene = chrome[i - 1];
			vector<WoCConnection>& currConnection = connectionsMap[gene.first];
			
			// Obtain the node and mask and iterate through said mask.
			Node* node = islands[gene.first];
			uint8 mask = gene.second;
			int bitCount = 0;
			for (mask; mask; mask >>= 1)
			{
				if (mask & 1)
				{
					// Obtain the direction and neighbor.
					Direction dir = static_cast<Direction>(bitCount % BITMASK_BOUNDARY);

					Node* nb = nullptr;
					vector<Neighbor*>::iterator nbiter = find_if(node->neighbors.begin(), node->neighbors.end(),
						[&dir](const Neighbor* n)
						{
							return n->neighborDirection == dir;
						});
					if (nbiter != node->neighbors.end())
					{
						nb = (*nbiter)->neighborNode;
					}

					if (nb)
					{
						// Obtain the actual connection.
						vector<WoCConnection>::iterator wociter = find_if(currConnection.begin(), currConnection.end(),
							[&nb](const WoCConnection& wocc)
							{
								return wocc.islandID == nb->nodeID;
							});
						if (wociter != currConnection.end())
						{
							// Increment if the connection is found.
							++(*wociter).count;
						}
						else
						{
							// Create a new one if not.
							currConnection.push_back({ nb->nodeID, bitCount < BITMASK_BOUNDARY, 1 });
						}
					}
				}
				++bitCount;
			}
		}
	}

	// Check to see the map filled.
	if (!connectionsMap.size())
	{
		return;
	}

	// Sort the map.
	for (pair<const int, vector<WoCConnection>>& connections : connectionsMap)
	{
		sort(connections.second.begin(), connections.second.end(),
			[](const WoCConnection& a, const WoCConnection& b)
			{
				return a.count > b.count;
			});
	}

	// Put the board back to its original state.
	ClearBridgesOnBoard();
	for (Node* node : islands)
	{
		UpdateBaseNeighborInfo(node);
	}

	// Create an empty chromosome.
	Chromosome emptyChrome;
	for (Node* node : islands)
	{
		emptyChrome.push_back(Gene(node->nodeID, 0));
	}

	// Make the wisdom chromosome.
	FitnessChromosome fWisdomChrome = FitnessChromosome(0.f, Chromosome(emptyChrome));
	Chromosome& wisdomChrome = fWisdomChrome.second;

	// Start building.
	for (Node* node : islands)
	{
		int id = node->nodeID;
		Gene& gene = wisdomChrome[id];
		int index = 0;
		// Iterate till the island is full or no connections are left.
		while(CalcConnectionsFromMask(gene.second) < node->value && index < connectionsMap[id].size())
		{
			WoCConnection wocc;
			if (!connectionsMap[id].empty())
			{
				bool bGeneVisited = true;
				// Same here.
				do
				{
					wocc = connectionsMap[id][index];
					if (find_if(wisdomChrome.begin(), wisdomChrome.end(),
						[&wocc](const Gene& g)
						{
							return g.first == wocc.islandID;
						}) != wisdomChrome.end())
					{
						bGeneVisited = false;
					}
				} while (++index < connectionsMap[id].size() && bGeneVisited);

				// Obtained the found connection neighbor.
				Neighbor* nb = nullptr;
				vector<Neighbor*>::iterator nbiter = find_if(node->neighbors.begin(), node->neighbors.end(),
					[&wocc](const Neighbor* n)
					{
						return n->neighborNode->nodeID == wocc.islandID;
					});
				if (nbiter != node->neighbors.end())
				{
					nb = *nbiter;
				}

				if (nb)
				{
					// Set the bits based off of direction.
					Direction dir = nb->neighborDirection;

					gene.second |= 1 << (static_cast<int>(dir) + (wocc.bSingleConnection ? 0 : BITMASK_BOUNDARY));
					wisdomChrome[nb->neighborNode->nodeID].second |= 1 << ((static_cast<int>(dir) ^ 1) + (wocc.bSingleConnection ? 0 : BITMASK_BOUNDARY));
				}
			}
		}
	}

	// If the chromosome is unique.
	if (CheckIfUnique(wisdomChrome))
	{
		// Fix, evaluate, and push.
		FixChromosomeConnections(wisdomChrome);
		EvaluateChromosome(fWisdomChrome);
		population.push_back(fWisdomChrome);
	}

	return;
}

bool HashiBoard::CheckIfUnique(const Chromosome& chromosome)
{
	for (const FitnessChromosome& fChrome : population)
	{
		int maxDupes = static_cast<int>(chromosome.size());
		const Chromosome& refChrome = fChrome.second;
		for (int i = 0; i < chromosome.size(); ++i)
		{
			const Gene& gene = chromosome[i];
			const Gene& refGene = refChrome[i];
			if (gene.second == refGene.second)
			{
				--maxDupes;
			}
		}

		if (!maxDupes)
		{
			return false;
		}
	}

	return true;
}

void HashiBoard::FixChromosomeConnections(Chromosome& chromosome)
{
	Chromosome shuffledChrome = chromosome;
	random_shuffle(shuffledChrome.begin(), shuffledChrome.end());
	// Fix mirroring connections between islands from different parents.
	FixMirroringConnections(shuffledChrome);
	// Fix excess connections that violate the island value.
	FixExcessConnections(shuffledChrome);

	for (Gene& gene : shuffledChrome)
	{
		Gene* origGene = nullptr;
		Chromosome::iterator iter = find_if(chromosome.begin(), chromosome.end(),
			[&gene](const Gene& g)
			{
				return g.first == gene.first;
			});
		if (iter != chromosome.end())
		{
			origGene = &(*iter);
		}

		if (origGene)
		{
			origGene->second = gene.second;
		}
	}
}

void HashiBoard::FixMirroringConnections(Chromosome& chromosome)
{
	for (Gene& gene : chromosome)
	{
		Node* node = islands[gene.first];
		//Before we start looking at mirroring connections, we need
		//to make sure that all the node's neighbors exist, even if
		//they might actually be hit by overlapping bridges.// Useful typedef for less typing.
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
		 
		//This will get checked and penalized later on, but for now, keep them.
		ClearBridgesOnBoard();
		UpdateBaseNeighborInfo(node);
		uint8& mask = gene.second;

		uint8 iterMask = mask;
		int bitCount = 0;
		// Iterate through the mask.
		for (int i = 0; i < 2 * BITMASK_BOUNDARY && iterMask; ++i)
		{
			if (iterMask & 1)
			{
				// Obtain the neighbor.
				Direction dir = static_cast<Direction>(bitCount % BITMASK_BOUNDARY);
				Neighbor* nb = nullptr;
				vector<Neighbor*>::iterator iter = find_if(node->neighbors.begin(), node->neighbors.end(),
					[&dir](const Neighbor* nb)
					{
						return nb->neighborDirection == dir;
					});
				if (iter != node->neighbors.end())
				{
					nb = *iter;
				}

				if (nb)
					// Obtain the gene.
				{
					int nbId = nb->neighborNode->nodeID;
					Gene* nbGene = nullptr;
					Chromosome::iterator citer = find_if(chromosome.begin(), chromosome.end(),
						[&nbId](const Gene& g)
						{
							return g.first == nbId;
						});
					if (citer != chromosome.end())
					{
						nbGene = &(*citer);
					}

					if (nbGene)
					{
						// Check to see if the mirrored bit is set. If not, set it but clear the single or double for the opposite connection.
						int numNbCon = (bitCount >= BITMASK_BOUNDARY ? 2 : 1);
						int shift = ((static_cast<int>(dir) ^ 1) + (numNbCon == 2 ? BITMASK_BOUNDARY : 0));
						if (!((nbGene->second >> shift) & 1))
						{
							if (nbGene->second >> (shift + (numNbCon == 2 ? -BITMASK_BOUNDARY : BITMASK_BOUNDARY)))
							{
								nbGene->second &= ~(1 << (shift + (numNbCon == 2 ? -BITMASK_BOUNDARY : BITMASK_BOUNDARY)));
							}
							nbGene->second |= 1 << shift;
						}
					}
				}
			}

			iterMask >>= 1;
			++bitCount;
		}
	}
}

void HashiBoard::FixExcessConnections(Chromosome& chromosome)
{
	for (Gene& gene : chromosome)
	{
		Node* node = islands[gene.first];
		uint8& mask = gene.second;

		// Ensure the number of connections does not exceed the node's value.
		while (CalcConnectionsFromMask(mask) > node->value)
		{
			// Randomly select a bit to turn off to reduce the number of connections.
			int bitToTurnOff = rand() % (2 * BITMASK_BOUNDARY);
			// Ensure the bit to turn off is set.
			while (!(mask & (1 << bitToTurnOff)) && bitToTurnOff < BITMASK_BOUNDARY)
			{
				bitToTurnOff = rand() % (2 * BITMASK_BOUNDARY);
			}
			mask &= ~(1 << bitToTurnOff);

			// Obtain the neighbor.
			Neighbor* nb = nullptr;
			vector<Neighbor*>::iterator iter = find_if(node->neighbors.begin(), node->neighbors.end(),
				[&bitToTurnOff](const Neighbor* nb)
				{
					Direction dir = static_cast<Direction>(bitToTurnOff % BITMASK_BOUNDARY);
					return nb->neighborDirection == dir;
				});
			if (iter != node->neighbors.end())
			{
				nb = *iter;
			}

			if (nb)
			{
				// Obtain the neighbor gene.
				int nbId = nb->neighborNode->nodeID;
				Gene* nbGene = nullptr;
				Chromosome::iterator citer = find_if(chromosome.begin(), chromosome.end(),
					[&nbId](const Gene& g)
					{
						return g.first == nbId;
					});
				if (citer != chromosome.end())
				{
					nbGene = &(*citer);
				}

				if (nbGene)
				{
					// Clear the bit.
					int shift = (bitToTurnOff < BITMASK_BOUNDARY ? ((bitToTurnOff % BITMASK_BOUNDARY) ^ 1) : (((bitToTurnOff % BITMASK_BOUNDARY) ^ 1) + BITMASK_BOUNDARY));
					nbGene->second &= ~(1 << shift);
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////
// For debuggin purposes for now, I'll delete this later. 
// Prints every single iteration to see the changes
////////////////////////////////////////////////////////////

void HashiBoard::PrintBoard() const {
	for (const auto& row : board) {
		for (const int cell : row) {
			if (cell > 0) {
				int bridgeCount = 0;
				for (int i = 0; i < BITMASK_BOUNDARY * 2; ++i) {
					if (cell & (1 << i)) {
						bridgeCount += (i < BITMASK_BOUNDARY) ? 2 : 1;
					}
				}
				cout << cell << " ";//<< "(" << bridgeCount << ") "; // Print island with bridge count to see where numbers are differing
			}
			else {
				// Print bridges
				switch (cell) {
				case HORI_SINGLE_BRIDGE: cout << "- "; break;
				case HORI_DOUBLE_BRIDGE: cout << "= "; break;
				case VERT_SINGLE_BRIDGE: cout << "| "; break;
				case VERT_DOUBLE_BRIDGE: cout << "||"; break;
				default: cout << ". ";
				}
			}
		}
		cout << endl;
	}
	cout << endl;
}

int HashiBoard::IsDisjoint() const {
	if (islands.empty()) return false;

	unordered_set<int> visited;
	stack<Node*> nodeStack;
	int disjointGroups = 0;

	for(Node* node : islands)
	{
		if (visited.find(node->nodeID) != visited.end()) continue;
		// Using DFS to check connectivity
		nodeStack.push(node);
		visited.insert(node->nodeID);

		while (!nodeStack.empty()) {
			Node* current = nodeStack.top();
			nodeStack.pop();

			for (Neighbor* neighbor : current->neighbors) {
				if (neighbor->numOfBridges > 0 && visited.find(neighbor->neighborNode->nodeID) == visited.end()) {
					visited.insert(neighbor->neighborNode->nodeID);
					nodeStack.push(neighbor->neighborNode);
				}
			}
		}
		++disjointGroups;
	}

	// Check if all nodes were visited
	return disjointGroups - 1;
}

void HashiBoard::ParseIntString(string intsString, vector<int>& ints)
{
	string tmp = intsString;
	stringstream ss(tmp);
	string token;
	while (getline(ss, token, ','))
	{
		int i = StringToInt(token);
		ints.push_back(i);
	}
}

void HashiBoard::ParseFloatString(string floatsString, vector<float>& floats)
{
	string tmp = floatsString;
	stringstream ss(tmp);
	string token;
	while (getline(ss, token, ','))
	{
		float f = StringToFloat(token);
		floats.push_back(f);
	}
}

int HashiBoard::StringToInt(string intString)
{
	int i = 0;
	int digits = 0;
	string str = intString;
	while (!str.empty())
	{
		char c = str.back(); str.pop_back();
		i += ASCII_ATOI(c) * (digits ? static_cast<int>(pow(10, digits)) : 1);
		++digits;
	}
	return i;
}

float HashiBoard::StringToFloat(string floatString)
{
	int whole = 0;
	int wholeCount = 0;
	int decimal = 0;
	int decimalCount = 0;
	string str = floatString;
	bool bDecimalFound = false;
	while (!str.empty())
	{
		char c = str.back(); str.pop_back();
		if (c == '.')
		{
			bDecimalFound = true;
			continue;
		}
		if (!bDecimalFound)
		{
			decimal += ASCII_ATOI(c) * (decimalCount ? static_cast<int>(pow(10, decimalCount)) : 1);
			++decimalCount;
		}
		else
		{
			whole += ASCII_ATOI(c) * (wholeCount ? static_cast<int>(pow(10, wholeCount)) : 1);
			++wholeCount;
		}
	}
	return static_cast<float>(whole) + (static_cast<float>(decimal) / static_cast<float>(pow(10, decimalCount)));
}

Direction HashiBoard::GetOppositeDirection(Direction currentDirection)
{
	switch (currentDirection)
	{
	case Direction::LEFT:	return Direction::RIGHT;
	case Direction::RIGHT:	return Direction::LEFT;
	case Direction::DOWN:	return Direction::UP;
	case Direction::UP:		return Direction::DOWN;
	default:				return Direction::INVALID;
	}
}
