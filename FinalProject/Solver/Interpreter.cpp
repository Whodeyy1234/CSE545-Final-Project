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
	SDL_DestroyTexture(texture);
	texture = SDL_CreateTexture(SDL->GetRenderer(), SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, SCREEN_WIDTH, SCREEN_HEIGHT);
	bLongerWidth = false;
	boardSizeX = 0;
	boardSizeY = 0;
	board.clear();
	islands.clear();

	return true;
}

bool HashiBoard::ParseBoardFile(string filePath)
{
	ifstream file(filePath, ios::in);
	if (file.is_open())
	{
		int numRows = 0;
		bool bBoardStarted = false;

		string line;
		getline(file, line);
		while (!file.eof())
		{
			getline(file, line);
			if (line.size())
			{
				const size_t openBracket = line.find('[');
				const size_t closeBracket = line.find(']');
				const size_t comma = line.find(',');
				if (openBracket != string::npos && !bBoardStarted)
				{
					bBoardStarted = true;
				}
				else if (openBracket != string::npos && bBoardStarted)
				{
					board.push_back(vector<int>());
				}
				else if (closeBracket != string::npos && bBoardStarted)
				{
					++numRows;
				}
				else if (comma != string::npos && bBoardStarted)
				{
					string::const_iterator sIter = line.cbegin();
					for (sIter; sIter != line.cend(); ++sIter)
					{
						int id = ASCII_ATOI(*sIter);
						if (id >= 0 && id <= 9)
						{
							board[numRows].push_back(id);
						}
					}
				}
			}
		}

		file.close();
	}

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
		for (index = col - 1; index > 0; index--) 
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
		for (index = row - 1; index > 0; index--)
		{
			parsedNodeValue = board[index][col];
			if (parsedNodeValue > 0) { return GetNodeAtCoords(index, col); }
			else if (parsedNodeValue < 0) { return nullptr; } //Bridge has been hit when a negative value. 
		}
		break;
	}
	return nullptr;
}

void HashiBoard::RenderBoard() 
{
	if (!board.size())
	{
		return;
	}

	// Swap render target to the member texture.
	SDL_SetRenderTarget(SDL->GetRenderer(), texture);
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