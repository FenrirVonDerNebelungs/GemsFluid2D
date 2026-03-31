// GemsFluid2DWin.cpp : Defines the entry point for the application.
//


#include "GemsFluid2DWin.h"
#include "GenImage.h"
//#include "Fluid2D.h"
#include "Test.h"

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
GenImage* g_pGenImage = nullptr; //class used to draw images
//Fluid2D* g_pFluid2D = nullptr; //class used to run the fluid simulation
Test* g_pTest = nullptr; //class used to run the fluid simulation

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
VOID                Release();
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_GEMSFLUID2DWIN, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_GEMSFLUID2DWIN));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    Release();
    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_GEMSFLUID2DWIN));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_GEMSFLUID2DWIN);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
	/*g_pFluid2D = new Fluid2D();
	s_WH grid_wh = g_pFluid2D->getGridWidthHeight();
	g_pFluid2D->launchCUDA();*/
	g_pTest = new Test();
	s_WH grid_wh = g_pTest->getGridWidthHeight();
    g_pTest->runTest();
    g_pGenImage = g_pTest->getTestImage();//new GenImage(grid_wh.width, grid_wh.height);
   //double* Ux = g_pTest->getUx();
   //g_pGenImage->genNormalizedImage(Ux);
   hInst = hInstance; // Store instance handle in our global variable
   s_WH blown_grid_wh = g_pTest->getBlownWidthHeight();
   int window_width = blown_grid_wh.width + 6;
   int window_height = blown_grid_wh.height + 6;
   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, window_width, window_height, nullptr, nullptr, hInstance, nullptr);//CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);


   return TRUE;
}
VOID Release()
{
    if (g_pTest != nullptr)
    {
        delete g_pTest;
        g_pTest = nullptr;
    }
    /*
    if(g_pFluid2D != nullptr)
    {
        delete g_pFluid2D;
        g_pFluid2D = nullptr;
	}*/
}
//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE: Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Parse the menu selections:
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
			BITMAPINFO* bmi = g_pGenImage->getBitmapInfo();
			s_WH wh = g_pGenImage->getWidthHeight();
			unsigned char* pImageData = g_pGenImage->getImageData();

            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            if (pImageData != nullptr) {
                SetDIBitsToDevice(hdc,
                    0, 0,               // xDest, yDest
                    wh.width, wh.height,     // nWidth, nHeight
                    0, 0,             // xSrc, ySrc
                    0,                // Start scan line
                    wh.height,          // number of scan lines
                    pImageData, // pointer to the array of color data
                    bmi,
                    DIB_RGB_COLORS);
            }
            EndPaint(hWnd, &ps);

        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_LBUTTONDOWN: {
        int x_click = LOWORD(lParam);
        int y_click = HIWORD(lParam);
        g_pGenImage = g_pTest->handleMouse();
        if(g_pTest->getMessage()!=nullptr)
            MessageBox(hWnd, g_pTest->getMessage(), L"Active", MB_OK);
		InvalidateRect(hWnd, nullptr, TRUE);/*invalidates the entire client area and causes a WM_PAINT message to be sent to the window procedure
                                              NULL means entire client area, true to erase background*/
        //MessageBox(hWnd, L"Left mouse button clicked", L"Mouse Click", MB_OK);
        break;
    }
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
