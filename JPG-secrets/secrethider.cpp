#include <iostream>
#include <atlimage.h>
#include <string>
using namespace std;

string charToBinary(char c);

int main()
{
	CImage img;
	cout << "Loading image...........";
	const char* sourcePath = "C:/Users/Steinar/source/repos/JPGEnveloper/x64/Debug/andean-condor-e1564404603879.jpg";
	img.Load(sourcePath);
	cout << "done" << endl;
	cout << "Loading secret..........";
	string secret = "Hello world!";
	secret = to_string(secret.length()) + ":" + secret;
	cout << "done" << endl;

	if (secret.length() > img.GetWidth() * img.GetHeight() * 8)
	{
		cout << "Image too small for message, must be at least " + to_string(img.GetWidth() * img.GetHeight() * 8) + " pixels" << endl;
	}
	else
	{
		cout << "Burying secret in jpg...";
		for (int i = 0; i < secret.length(); i++)
		{
			string binary = charToBinary(secret.at(i));
			for (int u = 0; u < 8; u++)
			{
				int pixelNo = i * 8 + u;
				int x = pixelNo % img.GetWidth();
				int y = pixelNo / img.GetWidth();
				COLORREF clr = img.GetPixel(x, y);

				char r = GetRValue(clr);
				char g = GetGValue(clr);
				char b = GetBValue(clr);

				if (b % 2 == 0 && binary[u] == '1') // partall - 0
				{
					b--;
					img.SetPixelRGB(x, y, r, g, b);
				}
				if (b % 2 == 1 && binary[u] == '1') // oddetall - 1
				{
					b++;
					img.SetPixelRGB(x, y, r, g, b);
				}
			}
		}
		cout << "done" << endl;
		cout << "Saving image............";
		img.Save("C:/Users/Steinar/source/repos/JPGEnveloper/x64/Debug/andean-condor-e1564404603879.jpg");
		cout << "done" << endl;
	}

	

	cout << "Task complete, press enter to exit";
	std::cin.ignore();
	return 0;
}

string charToBinary(char input)
{
	string output = "00000000";
	int diff = 128, c = input;

	for (int i = 0; i < 8; i++)
	{
		if (c >= diff)
		{
			output[i] = '1';
			c -= diff;
		}
		diff = diff / 2;
	}
	return output;
}