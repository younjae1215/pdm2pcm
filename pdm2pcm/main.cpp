#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define path "C:\\Users\\user\\Desktop\\pdm\\"

void PDM2PCM(std::vector<unsigned char>& pdmSignal, std::vector<int>& pcmSignal, int decimationFactor)
{
	int decimationF = decimationFactor / 8; // byte size
	int N = (pdmSignal.size() / decimationF); // PCM size, PCM size = (PDM size / decimation facor) ?		
	int k = 0;

	std::vector<int> buffer;

	/* Buffer PDM signal for further processing and decimation */
	for (int i = 0; i < N; i++)
	{
		int tmp = 0;
		for (int imod = 0; imod < decimationF; imod++)
		{
			unsigned char tmpBitSum = 0;
			unsigned char sig = pdmSignal.at(i * decimationF + imod);

			// endian check
			for (int iter = 0; iter < 8; ++iter)
			{				
				unsigned char tmpData = ((sig >> iter) & 1);
				tmpBitSum += tmpData;
			}
			tmp += tmpBitSum;
		}

		buffer.push_back(tmp);
	}

	pcmSignal.assign(N, 0);

	///* Low-Pass Filter Signal To Obtain PCM signal */
	//int filterParam = 10;
	//for (int i = 0; i < N - filterParam; i++)
	//{
	//	int tmpData = 0;
	//	for (int imod = 0; imod < filterParam; imod++)
	//	{
	//		tmpData += buffer[i + imod];
	//	}

	//	pcmSignal.at(i) = tmpData;
	//}

	int volume = 1;

	for (int i = 0; i < N - decimationFactor; i++)
	{
		int tmpData = 0;
		for (int imod = 0; imod < decimationFactor; imod++)
		{
			//pcmSignal[i] += pcmSignal[i + imod];
			tmpData += buffer[i + imod];
		}

		pcmSignal[i] = (tmpData) * volume;
	}
	
	return;
}

void write32(std::ofstream& os, int data)
{
	for (int i = 0; i < 4; ++i)
	{
		char tmpData = (data >> 8 * i);
		os << tmpData;
	}
}

void write16(std::ofstream& os, int data)
{
	for (int i = 0; i < 2; ++i)
	{
		char tmpData = (data >> 8 * i);
		os << tmpData;
	}
}

void makeWAV(int samplingRate, std::vector<int>& pcmSignal_left, std::vector<int>& pcmSignal_right, std::string& fileName)
{
	//https://blog.naver.com/PostView.nhn?blogId=psychoria&logNo=40139175382
	
	std::string tmpStr(path);
	tmpStr += "test.wav";

	std::ofstream writeFile;
	writeFile.open(tmpStr.c_str());
	writeFile << "RIFF";
	write32(writeFile, 36 + pcmSignal_left.size() * 2);
	writeFile << "WAVE";
	writeFile << "fmt ";
	write32(writeFile, 16);
	write16(writeFile, 1);
	write16(writeFile, 1 * 2);
	write32(writeFile, samplingRate);
	write32(writeFile, samplingRate * 2 * 2);
	write16(writeFile, 2 * 2);
	write16(writeFile, 16);
	writeFile << "data";
	write32(writeFile, pcmSignal_left.size() * 2);

	// write data
	for (int i = 0; i < pcmSignal_left.size(); ++i)
	{
		char data1 = pcmSignal_left.at(i);
		char data2 = pcmSignal_left.at(i) >> 8;
		char data3 = pcmSignal_right.at(i);
		char data4 = pcmSignal_right.at(i) >> 8;
		writeFile << data1;
		writeFile << data2;
		writeFile << data3;
		writeFile << data4;
	}

	writeFile.close();

	return;
}

void main()
{
	std::vector<unsigned char> pdmSignal_left;
	std::vector<unsigned char> pdmSignal_right;

	std::string fileName(path);
	fileName += "pdm.txt";

	std::ifstream readFile;
	readFile.open(fileName.c_str());

	if (readFile.is_open())
	{
		while (!readFile.eof())
		{
			std::string tmpStr;
			readFile >> tmpStr;
			int tmpData = std::atoi(tmpStr.c_str());

			unsigned char a = tmpData >> 8;
			unsigned char b = tmpData;

			pdmSignal_left.push_back(a);
			pdmSignal_right.push_back(b);
		}
	}
	readFile.close();

	std::vector<int> pcmSignal_left;
	std::vector<int> pcmSignal_right;

	int decimationFactor = 128;
	PDM2PCM(pdmSignal_left, pcmSignal_left, decimationFactor);
	PDM2PCM(pdmSignal_right, pcmSignal_right, decimationFactor);

	fileName.clear();
	fileName = path;
	fileName += "test.wav";

	int samplingRate = 16384;
	makeWAV(samplingRate, pcmSignal_left, pcmSignal_right, fileName);

	return;
}