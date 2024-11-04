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

	// ImGui Variables - Start ===============
	string filePath = "";
	string cachedFilePath = "";
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

		ImGui::Begin("Input");

		ImGui::InputText("File Path", &filePath);

		if (ImGui::Button("Show Board"))
		{
			if (cachedFilePath != filePath)
			{
				cachedFilePath = filePath;
				hashiBoard->Reset();
			}
			hashiBoard->Initialize(filePath);
		}

		ImGui::End();

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