#include <imgui.h>
#include <imgui_impl_sdl2.h>
#include <imgui_impl_sdlrenderer2.h>
#include <imgui_stdlib.h>
#include "Interpreter.h"

#include <SDL2Singleton.h>

SDLSingleton* SDLSingleton::mpInstance = nullptr;

HashiBoard* SetupHashiBoard() 
{
	HashiBoard* hashiBoard = new HashiBoard();
	hashiBoard->ParsePuzzle();
	return hashiBoard;
}

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

	HashiBoard* hashiBoard = SetupHashiBoard();

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

		ImGui::Begin("Hello, World!");

		ImGui::End();

		// Render
		ImGui::Render();

		SDL_RenderClear(SDL->GetRenderer());

		ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), SDL->GetRenderer());

		SDL_RenderPresent(SDL->GetRenderer());
	}

	// Shutdown ImGui
	ImGui_ImplSDLRenderer2_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();

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