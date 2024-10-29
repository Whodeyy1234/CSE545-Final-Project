#include "SDL2Singleton.h"

bool SDLSingleton::Initialize()
{
	// Initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
	{
		return false;
	}

	// Initialize SDL_ttf
	if (TTF_Init() == -1)
	{
		SDL_Quit();
		return false;
	}

	// Make the window
	mpWindow = SDL_CreateWindow(
		"TSP Solver",
		SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED,
		SCREEN_WIDTH,
		SCREEN_HEIGHT,
		SDL_WINDOW_SHOWN
	);

	// Quit if the window wasn't made
	if (!mpWindow)
	{
		TTF_Quit();
		SDL_Quit();
		return false;
	}

	// Make the renderer
	mpRenderer = SDL_CreateRenderer(mpWindow, -1, 0);

	// Quit if the renderer wasn't made
	if (!mpRenderer)
	{
		SDL_DestroyWindow(mpWindow);
		TTF_Quit();
		SDL_Quit();
		return false;
	}

	// Create a font
	mpFont = TTF_OpenFont(FONT_FILE_PATH, POINT_DIM * ID_SCALE);
	if (!mpFont)
	{
		SDL_DestroyRenderer(SDL->GetRenderer());
		SDL_DestroyWindow(SDL->GetWindow());
		TTF_Quit();
		SDL_Quit();
		return false;
	}

	// If everything initialized and made successfully, return true
	return true;
}

bool SDLSingleton::Destroy()
{
	// Destroys and quits everything SDL2 related.
	SDL_DestroyRenderer(SDL->GetRenderer());
	SDL_DestroyWindow(SDL->GetWindow());
	TTF_Quit();
	SDL_Quit();
	return true;
}