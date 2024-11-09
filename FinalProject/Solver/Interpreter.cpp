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
		for (index = col - 1; index >= 0; index--) 
		{
			parsedNodeValue = board[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::RIGHT:
		for (index = col + 1; index < boardSizeX; index++)
		{
			parsedNodeValue = board[row][index];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(row, index); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::DOWN:
		for (index = row + 1; index < boardSizeY; index++)
		{
			parsedNodeValue = board[index][col];
			if (parsedNodeValue > 0)  { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	case Direction::UP:
		for (index = row - 1; index >= 0; index--)
		{
			parsedNodeValue = board[index][col];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
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

bool HashiBoard::Update(Parameters params)
{
	// Process the algorithm.
	bool bFinished = Process(params);

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

	return bFinished;
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
		float perc = static_cast<float>(rand() / RAND_MAX);
		if (perc < params.crossoverProb)
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
		else if (perc < params.mutationProb)
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
		// @TODO: Put wisdom of crowds logic here.
	}
	// Outputting to csv if new best parent.
	if (population[0].first != bestPerc)
	{
		bestPerc = population[0].first;
		outputFile << currGen << "," << bestPerc << endl;
	}
	// Ending condition.
	if (bestPerc >= 1.0f || currGen > params.maxGenerations)
	{
		outputFile.close();
		return false;
	}

	return true;
}

bool HashiBoard::InitializePopulation(int populationSize, unsigned int seed)
{
	// Seed the random generator.
	cout << "Seed: " << seed << endl;
	srand(seed);

	outputFile.open("Output.csv");
	if (outputFile.is_open())
	{
		outputFile << "Seed," << seed << endl;
		outputFile << "Generation,Solved Perc" << endl;
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
		// Iterate through the genes and initialize connections.
		for (int i = 0; i < chromosome.size(); ++i)
		{
			Gene& gene = chromosome[i];
			InitializeIslandConnections(gene.first, gene.second, chromosome);
		}
		FixChromosomeConnections(chromosome);
		// Evaluate the fitness of the chromosome before leaving.
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
		if(nbGene)
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
					++numConnections;
					nbGene->second |= 1 << (static_cast<int>(nb->neighborDirection) ^ 1);
				}
			}
		}
	}
}

void HashiBoard::EvaluateChromosome(FitnessChromosome& chrome)
{
	// Fitness is defined as the completion of all island connections.
	int numRequiredConnections = 0;
	int numAcquiredConnections = 0;
	for (Gene& gene : chrome.second)
	{
		numRequiredConnections += islands[gene.first]->value;
		numAcquiredConnections += CalcConnectionsFromMask(gene.second);
	}

	chrome.first = static_cast<float>(numAcquiredConnections) / numRequiredConnections;
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
	SDL_Surface* TextSurface = TTF_RenderText_Solid(SDL->GetFont(), idString.data(), {255, 255, 255, 255});
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

void HashiBoard::PerformMutation(const vector<int>& mutationChromes) {
	for (int chromeIndex : mutationChromes) {
		Chromosome& chromosome = population[chromeIndex].second;

		// Randomly mutate genes in the chromosome
		for (Gene& gene : chromosome) {
			uint8& connection = gene.second;
			int bitToToggle = rand() % (2 * BITMASK_BOUNDARY);
			connection ^= (1 << bitToToggle);

			// Also toggle the corresponding bit in the neighbor's gene
			Node* node = islands[gene.first];
			for (Neighbor* neighbor : node->neighbors) {
				if ((bitToToggle % BITMASK_BOUNDARY == static_cast<int>(neighbor->neighborDirection)) ||
					(bitToToggle % BITMASK_BOUNDARY + BITMASK_BOUNDARY == static_cast<int>(neighbor->neighborDirection))) {
					Node* neighborNode = neighbor->neighborNode;
					for (Gene& neighborGene : chromosome) {
						if (neighborGene.first == neighborNode->nodeID) {
							neighborGene.second ^= (1 << ((bitToToggle + BITMASK_BOUNDARY) % (2 * BITMASK_BOUNDARY)));
							break;
						}
					}
				}
			}
		}

		FixChromosomeConnections(chromosome);
		EvaluateChromosome(population[chromeIndex]);
	}
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

void HashiBoard::ValidateConnections(Chromosome& chrome) {
	for (Gene& gene : chrome) {
		Node* node = islands[gene.first];
		int maxConnections = node->value; // max allowed connections for this node
		int currentConnections = CalcConnectionsFromMask(gene.second);

		if (currentConnections > maxConnections) {
			// Reduce excess connections by clearing bits in the mask
			ReduceExcessConnections(gene.second, currentConnections - maxConnections);
		}
	}
}

void HashiBoard::ReduceExcessConnections(uint8& mask, int excessConnections) {
	int bitCount = 0;
	while (excessConnections > 0 && bitCount < BITMASK_BOUNDARY * 2) {
		if (mask & (1 << bitCount)) { // Check if the bit is set
			mask &= ~(1 << bitCount); // Clear the bit
			excessConnections--;
		}
		bitCount++;
	}
}

void HashiBoard::FixChromosomeConnections(Chromosome& chromosome) 
{
	// Fix mirroring connections between islands from different parents.
	FixMirroringConnections(chromosome);
	// Fix overlapping connections.
	FixOverlappingConnections(chromosome);
	// Fix excess connections that violate the island value.
	FixExcessConnections(chromosome);
}

void HashiBoard::FixMirroringConnections(Chromosome& chromosome)
{
	for (Gene& gene : chromosome)
	{
		Node* node = islands[gene.first];
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

void HashiBoard::FixOverlappingConnections(Chromosome& chromosome)
{
	vector<pair<pair<int, int>, pair<int, int>>> ExistingConnections;
	for (const Gene& gene : chromosome)
	{
		Node* node = islands[gene.first];
		uint8 mask = gene.second;
		int bitCount = 0;
		for (mask; mask; mask >>= 1)
		{
			if (mask & 1)
			{
				Direction dir = static_cast<Direction>(bitCount % BITMASK_BOUNDARY);
				Node* nb = nullptr;
				vector<Neighbor*>::iterator iter = find_if(node->neighbors.begin(), node->neighbors.end(),
					[&dir](const Neighbor* nb)
					{
						return nb->neighborDirection == dir;
					});
				if (iter != node->neighbors.end())
				{
					nb = (*iter)->neighborNode;
				}

				if (nb)
				{
					pair<pair<int, int>, pair<int, int>> connection(
						pair<int, int>(node->coords[0], node->coords[1]), pair<int, int>(nb->coords[0], nb->coords[1])
					);

					if (find_if(ExistingConnections.begin(), ExistingConnections.end(),
						[&connection](const pair<pair<int, int>, pair<int, int>>& p)
						{
							return ((p.first.first == connection.first.first && p.first.second == connection.first.second &&
								p.second.first == connection.second.first && p.second.second == connection.second.second) ||
								(p.first.first == connection.second.first && p.first.second == connection.second.second &&
									p.second.first == connection.first.first && p.second.second == connection.first.second));
						}) == ExistingConnections.end())
					{
						ExistingConnections.push_back(connection);
					}
				}
			}
			++bitCount;
		}
	}

	for (Gene& gene : chromosome)
	{
		Node* node = islands[gene.first];
		uint8 mask = gene.second;
		int bitCount = 0;
		for (mask; mask; mask >>= 1)
		{
			if (mask & 1)
			{
				Direction dir = static_cast<Direction>(bitCount % BITMASK_BOUNDARY);
				Node* nb = nullptr;
				vector<Neighbor*>::iterator iter = find_if(node->neighbors.begin(), node->neighbors.end(),
					[&dir](const Neighbor* nb)
					{
						return nb->neighborDirection == dir;
					});
				if (iter != node->neighbors.end())
				{
					nb = (*iter)->neighborNode;
				}

				if (nb)
				{
					vector<pair<pair<int, int>, pair<int, int>>>::iterator eciter = find_if(ExistingConnections.begin(), ExistingConnections.end(),
						[&, dir, node, nb](const pair<pair<int, int>, pair<int, int>>& p)
						{
							Node* first = (dir == Direction::LEFT || dir == Direction::UP) ? nb : node;
							Node* second = (first == node) ? nb : node;

							int width = second->coords[1] - first->coords[1];
							int height = second->coords[0] - first->coords[0];

							pair<int, int> cFirst = (p.first.first <= p.second.first && p.first.second <= p.second.second) ? p.first : p.second;
							pair<int, int> cSecond = (cFirst == p.first) ? p.second : p.first;

							int cWidth = cSecond.second - cFirst.second;
							int cHeight = cSecond.first - cFirst.first;

							if (width == cWidth && height == cHeight) return false;

							if (!width) ++width;
							if (!height) ++height;

							if (!cWidth) ++cWidth;
							if (!cHeight) ++cHeight;

							return (first->coords[1] < cFirst.second + cWidth &&
								first->coords[1] + width > cFirst.second &&
								first->coords[0] < cFirst.first + cHeight &&
								first->coords[0] + height > cFirst.first);
						});
					if (eciter != ExistingConnections.end())
					{
						Gene* nbGene = nullptr;
						Chromosome::iterator giter = find_if(chromosome.begin(), chromosome.end(),
							[&nb](const Gene& g)
							{
								return g.first == nb->nodeID;
							});
						if (giter != chromosome.end())
						{
							nbGene = &(*giter);
						}

						if(nbGene)
						{
							gene.second &= ~(1 << bitCount);
							nbGene->second &= ~(1 << ((static_cast<int>(dir) ^ 1) + (bitCount >= BITMASK_BOUNDARY ? BITMASK_BOUNDARY : 0)));
							vector<pair<pair<int, int>, pair<int, int>>>::iterator citer = find(ExistingConnections.begin(), ExistingConnections.end(),
								pair<pair<int, int>, pair<int, int>>(pair<int, int>(node->coords[0], node->coords[1]), pair<int, int>(nb->coords[0], nb->coords[1])));
							if (citer != ExistingConnections.end())
							{
								ExistingConnections.erase(citer);
							}
						}
					}
				}
			}
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
			while (!(mask & (1 << bitToTurnOff)))
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
					nbGene->second &= ~(1 << (bitToTurnOff < BITMASK_BOUNDARY ? (bitToTurnOff % BITMASK_BOUNDARY) ^ 1 : (bitToTurnOff % BITMASK_BOUNDARY) ^ 1 << BITMASK_BOUNDARY));
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
				cout << cell << "(" << bridgeCount << ") "; // Print island with bridge count to see where numbers are differing
			}
			else {
				// Print bridges
				switch (cell) {
				case HORI_SINGLE_BRIDGE: cout << "- "; break;
				case HORI_DOUBLE_BRIDGE: cout << "= "; break;
				case VERT_SINGLE_BRIDGE: cout << "| "; break;
				case VERT_DOUBLE_BRIDGE: cout << "|| "; break;
				default: cout << ". ";
				}
			}
		}
		cout << endl;
	}
	cout << endl;
}
