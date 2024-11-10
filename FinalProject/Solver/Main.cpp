#include <imgui.h>
#include <imgui_impl_sdl2.h>
#include <imgui_impl_sdlrenderer2.h>
#include <imgui_stdlib.h>
#include "Interpreter.h"

#include <SDL2Singleton.h>

SDLSingleton* SDLSingleton::mpInstance = nullptr;

int main(int argc, char* argv[])
{
	// Initialize all things SDL.
	if (!SDL->Initialize())
	{
		return -1;
	}

	// ImGui initialization.
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui_ImplSDL2_InitForSDLRenderer(
		SDL->GetWindow(),
		SDL->GetRenderer()
	);
	ImGui_ImplSDLRenderer2_Init(SDL->GetRenderer());

	HashiBoard* hashiBoard = new HashiBoard();
	Parameters params;
	BatchParams bparams;

	// ImGui Variables - Start ===============
	string filePath = "";					// File path to the board.
	// Algorithm Parameters - Start ==========
	bool bWithSeed = false;					// Whether or not to seed the random generator.
	int seed = static_cast<int>(time(0));	// Seed of the random generator.
	int populationSize = 200;				// Size of the population.
	float crossoverProb = 0.7f;				// Probability a crossover occurs.
	float mutationProb = 0.05f;				// Probability a mutation occurs.
	int maxGenerations = 2000;				// The max generations the algorithm runs.
	bool bWithWisdom = false;				// Whether or not to use Wisdom of Crowds
	int gensPerWisdom = 50;					// The number of generations between each wisdom path.
	float elitismPerc = 0.3f;				// The percentage of best paths to use.
	// Algorithm Parameters - End ============
	// Batch Testing Parameters - Start ======
	bool bWithBatch = false;
	int numExectuions = 1;
	string populationSizes = "";
	string crossoverProbs = "";
	string mutationProbs = "";
	string smaxGenerations = "";
	string gensPerWisdoms = "";
	string elitismPercs = "";
	// Batch Testing Parameters - End ========
	bool bAlgRunning = false;
	// ImGui Variables - End =================

	// Main loop.
	bool bRunning = true;
	while (bRunning)
	{
		// Main input loop.
		SDL_Event e;
		while (SDL_PollEvent(&e))
		{
			if (e.type == SDL_QUIT)
			{
				bRunning = false;
			}
			ImGui_ImplSDL2_ProcessEvent(&e);
		}

		// New ImGui frame.
		ImGui_ImplSDLRenderer2_NewFrame();
		ImGui_ImplSDL2_NewFrame();
		ImGui::NewFrame();

		if (!bAlgRunning)
		{
			ImGui::Begin("Input");

			ImGui::InputText("File Path", &filePath);

			if (ImGui::Button("Show Board"))
			{
				hashiBoard->Reset();
				hashiBoard->Initialize(filePath);
			}

			if (ImGui::CollapsingHeader("Genetic Parameters", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::Checkbox("Set Seed", &bWithSeed);
				if (bWithSeed)
				{
					ImGui::InputInt("Seed", &seed);
				}
				ImGui::InputInt("Population Size", &populationSize);
				ImGui::InputFloat("Crossover Probability", &crossoverProb);
				ImGui::InputFloat("Mutation Probability", &mutationProb);
				ImGui::InputInt("Maximum Generations", &maxGenerations);
				ImGui::Checkbox("With Wisdom", &bWithWisdom);
				if (bWithWisdom)
				{
					ImGui::InputInt("Generations Per Wisdom", &gensPerWisdom);
					ImGui::InputFloat("Elitism Percentage", &elitismPerc);
				}
			}

			if (ImGui::CollapsingHeader("Batch Testing Parameters", ImGuiTreeNodeFlags_None))
			{
				ImGui::Checkbox("Perform Batch Testing", &bWithBatch);
				if (bWithBatch)
				{
					ImGui::Text("Input Format: #,#,#,...");
					ImGui::InputText("Population Sizes", &populationSizes);
					ImGui::InputText("Crossover Probabilities", &crossoverProbs);
					ImGui::InputText("Mutation Probabilities", &mutationProbs);
					ImGui::InputText("Maximum Generations", &smaxGenerations);
					ImGui::InputText("Generations Per Wisdoms", &gensPerWisdoms);
					ImGui::InputText("Elitism Percentages", &elitismPercs);
				}
			}

			if (ImGui::Button("Run"))
			{
				if (!bWithBatch)
				{
					params = { static_cast<unsigned int>(bWithSeed ? seed : time(0)), populationSize, crossoverProb, mutationProb, maxGenerations, bWithWisdom, gensPerWisdom, elitismPerc };
				}
				else
				{
					bparams.puzzleFilePath = filePath;
					HashiBoard::ParseIntString(populationSizes, bparams.populationSizes);
					HashiBoard::ParseFloatString(crossoverProbs, bparams.crossoverProbs);
					HashiBoard::ParseFloatString(mutationProbs, bparams.mutationProbs);
					HashiBoard::ParseIntString(smaxGenerations, bparams.maxGenerations);
					HashiBoard::ParseIntString(gensPerWisdoms, bparams.gensPerWisdoms);
					HashiBoard::ParseFloatString(elitismPercs, bparams.elitismPercs);
				}
				bAlgRunning = true;
			}

			ImGui::End();
		}
		else
		{
			ImGui::Begin("Quitter");
			if (ImGui::Button("Quit"))
			{
				bAlgRunning = false;
			}
			else
			{
				bAlgRunning = bWithBatch ? hashiBoard->BatchUpdate(bparams) : hashiBoard->Update(params);
			}
			ImGui::End();
		}

		// Render
		ImGui::Render();

		SDL_RenderClear(SDL->GetRenderer());

		if (hashiBoard)
		{
			hashiBoard->RenderBoard();
		}

		ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), SDL->GetRenderer());

		SDL_RenderPresent(SDL->GetRenderer());
	}

	// Shutdown ImGui
	ImGui_ImplSDLRenderer2_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();

	delete hashiBoard;

	// Shutdown SDL.
	if (SDL->Destroy())
	{
		delete SDL;
		return 0;
	}
	else
	{
		delete SDL;
		return -1;
	}
}