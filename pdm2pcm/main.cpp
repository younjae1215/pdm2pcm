#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "OpenPDMFilter.h"

#define path "C:\\Users\\user\\Desktop\\pdm\\"
#define MAX_DECIMATION_FACTOR 128
#define useOpenPDMLib 1

// 256, 32, 32
void pdmFilter(uint16_t *PDMBuf, uint16_t *PCMBuf, uint32_t pcm_buff_size)
{
	uint32_t i = 0, j = 0, k = 0;
	uint32_t Current_Value = 0, Prev_Value = 0;
	uint32_t pcm_data = 0;
	uint8_t num_samples = 16;
	uint8_t decimationFactor = MAX_DECIMATION_FACTOR;
	uint8_t decimationF = decimationFactor / 16; // word size
	static uint32_t prev_buffer[8];


	for (k = 0; k < pcm_buff_size; k++)	// 0~31
	{
		pcm_data = 0;

		/* down sampling */
		for (j = 0; j < decimationF; j++)	// 0~7
		{
			Current_Value = PDMBuf[j + (k * decimationF)];
			Prev_Value = prev_buffer[j];

			/*Copy from the Current Data Sample to Previous Data Sample*/
			prev_buffer[j] = Current_Value;
			for (i = 0; i < num_samples; i++)	// 0~15
			{
				pcm_data += ((((Prev_Value & 0x8000) && 1) * (i + (num_samples * j))) + ((Current_Value & 0x8000) && 1) * (MAX_DECIMATION_FACTOR - i - (num_samples * j)));
				Current_Value = Current_Value << 1;
				Prev_Value = Prev_Value << 1;
			}
		}
		/* Now Perform the Shift - Since already in the 4096 range for Order 7 shift only by 4 and so on*/
//		PCMBuf[k] = (int16_t)( (pcm_data << 4) - 32768);//pcm_data;//(int16_t)( (pcm_data << 4) - 32768);
		PCMBuf[k] = (int16_t)((pcm_data) << 4);
		//		PCMBuf[k] = (int16_t)(pcm_data);
	}
}

void PDM2PCM(std::vector<unsigned char>& pdmSignal, std::vector<uint16_t>& pcmSignal, int decimationFactor)
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

	/* Low-Pass Filter Signal To Obtain PCM signal */
	int filterParam = 100;
	for (int i = 0; i < N - filterParam; i++)
	{
		int tmpData = 0;
		for (int imod = 0; imod < filterParam; imod++)
		{
			tmpData += buffer[i + imod];
		}

		buffer.at(i) = tmpData;
	}

	float volume = 1;

	for (int i = 0; i < N - decimationFactor; i++)
	{
		int tmpData = 0;
		for (int imod = 0; imod < decimationFactor; imod++)
		{
			//pcmSignal[i] += pcmSignal[i + imod];
			tmpData += buffer[i + imod];
		}

		pcmSignal[i] = tmpData * volume;
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

void makeWAV(int samplingRate, std::vector<int16_t>& pcmSignal_left, std::vector<int16_t>& pcmSignal_right)
{
	//https://blog.naver.com/PostView.nhn?blogId=psychoria&logNo=40139175382
	
	std::string tmpStr(path);
	tmpStr += "test.wav";

	std::ofstream writeFile;
	writeFile.open(tmpStr.c_str());
	writeFile << "RIFF";
	write32(writeFile, 36 + pcmSignal_left.size() * 2 * 2);
	writeFile << "WAVE";
	writeFile << "fmt ";
	write32(writeFile, 16);
	write16(writeFile, 1);
	write16(writeFile, 1 * 2);
	write32(writeFile, samplingRate);			// sample rate, 16K
	write32(writeFile, samplingRate * 2 * 2);	// byte rate, 16K * 2 channel * (16 / 8)
	write16(writeFile, 2 * 2);
	write16(writeFile, 16);		// bitPerSample, 16bit
	writeFile << "data";
	write32(writeFile, pcmSignal_left.size() * 2 * 2);	// subchunk2 size, (bitPerSample / 8) * 2 channel * sample numbers

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

void makeWAV_uint16(int samplingRate, std::vector<uint16_t>& pcmSignal, int channel)
{
	//https://blog.naver.com/PostView.nhn?blogId=psychoria&logNo=40139175382

	std::string tmpStr(path);
	tmpStr += "test_.wav";

	std::ofstream writeFile;
	writeFile.open(tmpStr.c_str());
	writeFile << "RIFF";
	write32(writeFile, 36 + pcmSignal.size());
	writeFile << "WAVE";
	writeFile << "fmt ";
	write32(writeFile, 16);
	write16(writeFile, 1);
	write16(writeFile, 1 * channel);					// channel
	write32(writeFile, samplingRate);			// sample rate, 16K
	write32(writeFile, samplingRate * 2 * channel);	// byte rate, 16K * 2 channel * (16 / 8)
	write16(writeFile, 2 * channel);
	write16(writeFile, 16);		// bitPerSample, 16bit
	writeFile << "data";
	write32(writeFile, pcmSignal.size());	// subchunk2 size, (bitPerSample / 8) * 2 channel * sample numbers

	// write data
	for (int i = 0; i < pcmSignal.size(); ++i)
	{
		write16(writeFile, pcmSignal.at(i));
	}

	writeFile.close();

	return;
}

void main()
{
	std::string fileName(path);
	//fileName += "200916_sample\\BlueCoin_Log_PDM_N002.txt";
	fileName += "BlueCoin_Log_PDM_N002.txt";

	if (!useOpenPDMLib)
	{
		std::ifstream readFile;
		readFile.open(fileName.c_str());

		std::vector<unsigned char> pdmSignal;

		if (readFile.is_open())
		{
			int i = 7;
			uint8_t tmpValue = 0;
			while (!readFile.eof())
			{
				std::string tmpStr;
				readFile >> tmpStr;
				uint16_t tmpData = std::atoi(tmpStr.c_str());

				tmpValue += tmpData * pow(2, i % 8);

				//if (i != 0 && i % 8 == 0)
				if (i == 0)
				{
					pdmSignal.push_back(tmpValue);
					tmpValue = 0;
					i = 7;
					continue;
				}
				--i;
			}
		}
		readFile.close();

		//std::vector<char> pdmSignal;

		//if (readFile.is_open())
		//{
		//	readFile.seekg(0, std::ios::end);
		//	int sz = readFile.tellg();
		//	readFile.seekg(0, std::ios::beg);

		//	pdmSignal.assign(sz, 0);
		//	readFile.read(pdmSignal.data(), sz);
		//}

		//std::vector<uint16_t> pdmSignal_uint16;

		//if (readFile.is_open())
		//{
		//	while (!readFile.eof())
		//	{
		//		std::string tmpStr;
		//		readFile >> tmpStr;
		//		uint16_t tmpData = std::atoi(tmpStr.c_str());

		//		unsigned char a = tmpData >> 8;
		//		unsigned char b = tmpData;

		//		pdmSignal_left.push_back(a);
		//		pdmSignal_right.push_back(b);

		//		pdmSignal_uint16.push_back(tmpData);
		//	}
		//}

		//readFile.close();

		//std::vector<unsigned char> pdmSignal_left;
		//std::vector<unsigned char> pdmSignal_right;

		//for (int i = 0; i < pdmSignal.size(); ++i)
		//{
		//	if (i % 2 == 0)
		//		pdmSignal_left.push_back(pdmSignal.at(i));
		//	else
		//		pdmSignal_right.push_back(pdmSignal.at(i));
		//}

		//for (int i = 0; i < pdmSignal_left.size(); ++i)
		//{
		//	uint16_t tmpData = (pdmSignal_left.at(i) << 8) + pdmSignal_right.at(i);
		//	pdmSignal_uint16.push_back(tmpData);
		//}

		std::vector<uint16_t> pcmSignal_left;
		//std::vector<int16_t> pcmSignal_right;

		int decimationfactor = 128;
		//PDM2PCM(pdmSignal_left, pcmSignal_left, decimationfactor);
		//PDM2PCM(pdmSignal_right, pcmSignal_right, decimationfactor);
		
		PDM2PCM(pdmSignal, pcmSignal_left, decimationfactor);

		//std::vector<uint16_t> pcmsignal;
		//uint16_t pcm_buff_size = 32;

		//int size = pdmSignal.size() / 256;

		//for (int i = 0; i < size - 1; ++i)
		//{
		//	std::vector<uint16_t> pdmbuff(pdmSignal.begin() + i * 256, pdmSignal.begin() + ((i + 1) * 256));
		//	std::vector<uint16_t> pcmbuff(pcm_buff_size, 0);
		//	pdmFilter(pdmbuff.data(), pcmbuff.data(), pcm_buff_size);

		//	pcmsignal.insert(pcmsignal.end(), pcmbuff.begin(), pcmbuff.end());
		//}

		//int samplingrate = 16000;
		int samplingrate = 8000;
		//makeWAV(samplingrate, pcmSignal_left, pcmSignal_right);
		makeWAV_uint16(samplingrate, pcmSignal_left, 1);
	}
	else
	{
		std::ifstream readFile;
		readFile.open(fileName.c_str());

		std::vector<char> pdmSignal;

		//if (readFile.is_open())
		//{
		//	int i = 7;
		//	uint8_t tmpValue = 0;
		//	while (!readFile.eof())
		//	{
		//		std::string tmpStr;
		//		readFile >> tmpStr;
		//		uint16_t tmpData = std::atoi(tmpStr.c_str());

		//		tmpValue += tmpData * pow(2, i % 8);

		//		//if (i != 0 && i % 8 == 0)
		//		if (i == 0)
		//		{
		//			pdmSignal.push_back(tmpValue);
		//			tmpValue = 0;
		//			i = 7;
		//			continue;
		//		}
		//		--i;
		//	}
		//}
		//readFile.close();

		if (readFile.is_open())
		{
			readFile.seekg(0, std::ios::end);
			int sz = readFile.tellg();
			readFile.seekg(0, std::ios::beg);

			pdmSignal.assign(sz, 0);
			readFile.read(pdmSignal.data(), sz);
		}

		unsigned int pdmSamplingF, decimationF, pcmSamplingF, pdmBufLen, pcmBufLen;
		TPDMFilter_InitStruct filter;

		pdmSamplingF = 1024000 * 2;
		decimationF = 128;
		pcmSamplingF = pdmSamplingF / decimationF;

		pdmBufLen = pdmSamplingF / 1000;		// 1024
		pcmBufLen = pdmBufLen / decimationF;	// 8

		/* Initialize Open PDM library */
		filter.Fs = pcmSamplingF;
		filter.nSamples = pcmBufLen;
		filter.LP_HZ = pcmSamplingF / 2;
		filter.HP_HZ = 10;
		filter.In_MicChannels = 1;
		filter.Out_MicChannels = 1;
		filter.Decimation = decimationF;
		Open_PDM_Filter_Init(&filter);

		int finished = 0;
		int dataCount = 0;
		int ret = 0;

		std::vector<uint16_t> pcmSignal;
		int pdmBuffSize = pdmBufLen / 8;
		int pcmBuffSize = pdmBufLen / decimationF;
		int size = pdmBuffSize;

		for (int i = 0; i < pdmSignal.size() / size; ++i)
		{
			std::vector<uint8_t> pdmBuff(pdmSignal.begin() + i * size, pdmSignal.begin() + ((i + 1) * size));
			std::vector<int16_t> pcmBuff(pcmBuffSize, 0);
			Open_PDM_Filter_128(pdmBuff.data(), pcmBuff.data(), 1, &filter);

			pcmSignal.insert(pcmSignal.end(), pcmBuff.begin(), pcmBuff.end());
		}

		makeWAV_uint16(pcmSamplingF, pcmSignal, 2);
	}

	return;
}