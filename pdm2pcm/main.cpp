#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "OpenPDMFilter.h"

#define path "C:\\Users\\user\\Desktop\\pdm\\"
#define MAX_DECIMATION_FACTOR 128

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

	int volume = 10;

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

void makeWAV(int samplingRate, std::vector<int>& pcmSignal_left, std::vector<int>& pcmSignal_right)
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

void makeWAV_uint16(int samplingRate, std::vector<int16_t>& pcmSignal)
{
	//https://blog.naver.com/PostView.nhn?blogId=psychoria&logNo=40139175382

	std::string tmpStr(path);
	tmpStr += "test_.wav";

	std::ofstream writeFile;
	writeFile.open(tmpStr.c_str());
	writeFile << "RIFF";
	write32(writeFile, 36 + pcmSignal.size() * 2);
	writeFile << "WAVE";
	writeFile << "fmt ";
	write32(writeFile, 16);
	write16(writeFile, 1);
	write16(writeFile, 1 * 2);					// channel
	write32(writeFile, samplingRate);			// sample rate, 16K
	write32(writeFile, samplingRate * 2 * 2);	// byte rate, 16K * 2 channel * (16 / 8)
	write16(writeFile, 2 * 2);
	write16(writeFile, 16);		// bitPerSample, 16bit
	writeFile << "data";
	write32(writeFile, pcmSignal.size() * 2);	// subchunk2 size, (bitPerSample / 8) * 2 channel * sample numbers

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
	std::vector<unsigned char> pdmSignal_left;
	std::vector<unsigned char> pdmSignal_right;

	std::vector<uint8_t> pdmSignal;

	std::string fileName(path);
	fileName += "bellazio.txt";

	std::ifstream readFile;
	readFile.open(fileName.c_str());

	if (readFile.is_open())
	{
		int i = 0;
		uint8_t tmpValue = 0;
		while (!readFile.eof())
		{
			std::string tmpStr;
			readFile >> tmpStr;
			uint16_t tmpData = std::atoi(tmpStr.c_str());

			tmpValue += tmpData * pow(2, i % 8);

			if (i != 0 && i % 8 == 0)
			{
				unsigned char a = tmpData >> 8;
				unsigned char b = tmpData;

				pdmSignal_left.push_back(a);
				pdmSignal_right.push_back(b);

				//pdmSignal.push_back(tmpData);
				pdmSignal.push_back(tmpValue);
				tmpValue = 0;
			}
			++i;
		}
	}
	readFile.close();

	//std::vector<int> pcmSignal_left;
	//std::vector<int> pcmSignal_right;

	//int decimationFactor = 128;
	//PDM2PCM(pdmSignal_left, pcmSignal_left, decimationFactor);
	//PDM2PCM(pdmSignal_right, pcmSignal_right, decimationFactor);
	//
	//std::vector<uint16_t> pcmSignal;
	//uint16_t pcm_buff_size = 32;
	//int size = pdmSignal.size() / 256;

	//for (int i = 0; i < size - 1; ++i)
	//{
	//	std::vector<uint16_t> pdmBuff(pdmSignal.begin() + i * 256, pdmSignal.begin() + ((i + 1) * 256));
	//	std::vector<uint16_t> pcmBuff(pcm_buff_size, 0);
	//	pdmFilter(pdmBuff.data(), pcmBuff.data(), pcm_buff_size);

	//	pcmSignal.insert(pcmSignal.end(), pcmBuff.begin(), pcmBuff.end());
	//}

	//int samplingRate = 8000;
	//makeWAV(samplingRate, pcmSignal_left, pcmSignal_right);
	//makeWAV_uint16(samplingRate, pcmSignal);

	unsigned int pdmSamplingF, decimationF, pcmSamplingF, pdmBufLen, pcmBufLen;
	//uint8_t* pdmBuf;
	//int16_t* pcmBuf;
	TPDMFilter_InitStruct filter;

	pdmSamplingF = 1024000;
	decimationF = 128;
	pcmSamplingF = pdmSamplingF / decimationF;

	pdmBufLen = pdmSamplingF / 1000;		// 1024
	pcmBufLen = pdmBufLen / decimationF;	// 8

	//pdmBuf = malloc(pdmBufLen / 8);
	//std::vector<uint8_t> pdmBuf(pdmBufLen / 8, 0);

	//pcmBuf = malloc(sizeof(int16_t)*pcmBufLen);
	//std::vector<uint8_t> pdmBuf(pcmBufLen, 0);

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

	//while (finished == 0) {
	//	/* Grab 1ms data from stdin */
	//	dataCount = 0;
	//	while ((dataCount < pdmBufLen / 8) && (finished == 0)) {
	//		ret = read(STDIN_FILENO, pdmBuf + dataCount, pdmBufLen / 8 - dataCount);
	//		if (ret < 0) {
	//			fprintf(stderr, "Error reading from STDIN: %s\n", strerror(errno));
	//			exit(errno);
	//		}

	//		if (ret == 0) {
	//			fprintf(stderr, "Decoding complete!\n");
	//			finished = 1;
	//		}

	//		dataCount += ret;
	//	}

	//	/* Decode PDM. Oldest PDM bit is MSB */
	//	Open_PDM_Filter_128(pdmBuf, pcmBuf, 1, &filter);

	//	/* Emit PCM decoded data to stdout */
	//	dataCount = 0;
	//	while (dataCount < sizeof(int16_t)*pcmBufLen) {
	//		ret = write(STDOUT_FILENO, pcmBuf + dataCount, sizeof(int16_t)*pcmBufLen - dataCount);
	//		if (ret < 0) {
	//			fprintf(stderr, "Error writing to STDOUT: %s\n", strerror(errno));
	//			exit(errno);
	//		}

	//		dataCount += ret;
	//	}
	//}

	std::vector<int16_t> pcmSignal;
	int pdmBuffSize = pdmBufLen / 8;
	int pcmBuffSize = pdmBufLen / decimationF;
	int size = 8;

	for (int i = 0; i < pdmSignal.size() / size; ++i)
	{
		std::vector<uint8_t> pdmBuff(pdmSignal.begin() + i * size, pdmSignal.begin() + ((i + 1) * size));
		std::vector<int16_t> pcmBuff(size, 0);
		Open_PDM_Filter_128(pdmBuff.data(), pcmBuff.data(), 1, &filter);

		pcmSignal.insert(pcmSignal.end(), pcmBuff.begin(), pcmBuff.end());
	}

	makeWAV_uint16(pcmSamplingF, pcmSignal);

	return;
}