#pragma once

#include <SDL.h>
#include <SDL_ttf.h>

#define SDL SDLSingleton::GetInstance()
#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define FONT_FILE_PATH "../Shared/Lib/SDL2/VCR_OSD_MONO_1.001.ttf"
#define POINT_DIM 5
#define ID_SCALE 4

/// <summary>
/// SDL2 Global Singleton for obtaining SDL2 related attributes.
/// </summary>
class SDLSingleton
{
public:
	/// <summary>
	/// Instance accessor.
	/// </summary>
	/// <returns>An instance of the SDLSingleton.</returns>
	static SDLSingleton* GetInstance()
	{
		if (!mpInstance)
		{
			mpInstance = new SDLSingleton();
		}

		return mpInstance;
	}

	/* Construction and Destruction */
	/// <summary>
	/// Initalizes the SDLSingleton.
	/// </summary>
	/// <returns>Successful initialization or not.</returns>
	bool Initialize();
	/// <summary>
	/// Destroys the SDLSingleton.
	/// </summary>
	/// <returns>Successful destruction or not.</returns>
	bool Destroy();

	/* Accessors */
	/// <summary>
	/// Obtains the SDL2 window.
	/// </summary>
	/// <returns>The SDL2 window.</returns>
	SDL_Window* GetWindow() { return mpWindow; }
	/// <summary>
	/// Obtains the SDL2 renderer.
	/// </summary>
	/// <returns>The SDL2 renderer.</returns>
	SDL_Renderer* GetRenderer() { return mpRenderer; }
	/// <summary>
	/// Obtains the SDL2 font.
	/// </summary>
	/// <returns>The SDL2 font.</returns>
	TTF_Font* GetFont() { return mpFont; }

	/* Private Methods */
private:
	/// <summary>
	/// Private Constructor
	/// </summary>
	SDLSingleton() : mpWindow(nullptr), mpRenderer(nullptr), mpFont(nullptr) {};

	/* Attributes */
private:
	/// <summary>
	/// The SDLSingleton instance.
	/// </summary>
	static SDLSingleton* mpInstance;
	/// <summary>
	/// The window instance.
	/// </summary>
	SDL_Window* mpWindow;
	/// <summary>
	/// The renderer instance.
	/// </summary>
	SDL_Renderer* mpRenderer;
	/// <summary>
	/// The font instance.
	/// </summary>
	TTF_Font* mpFont;
};