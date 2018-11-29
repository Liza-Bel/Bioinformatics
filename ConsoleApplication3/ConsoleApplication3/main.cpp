#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

vector<int> sort(vector<int> s) {
	int tmp = 0;
	int j;
	for (int i = 0; i < s.size(); i++) {
		j = i;
		for (int k = i; k < s.size(); k++) {
			if (s[j] > s[k]) {
				j = k;
			}
		}
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	}
	return s;
}

vector<int> generating_theoretical_spectrum(string pepide) {
	int n = 0;
	int size = pepide.length()*(pepide.length() - 1) + 2;
	vector<int> arr;
	arr.resize(size);

	for (int i = 0; i < size; i++) {
		arr[i] = 0;
	}
	int k = 1;
	int sum = 0;
	string tmp;
	for (int j = 1; j < pepide.length() + 1; j++) {
		for (int m = 0; m < pepide.length(); m++) {

			if (j == pepide.length())
				tmp = pepide;

			else if (m + j > pepide.length()) {
				tmp = pepide.substr(m, pepide.length() - m);
				tmp += pepide.substr(0, j + m - pepide.length());
			}
			else tmp = pepide.substr(m, j);

			for (int i = 0; i < tmp.length(); i++) {
				if (tmp[i] == 'G')
					n = 57;
				if (tmp[i] == 'A')
					n = 71;
				if (tmp[i] == 'S')
					n = 87;
				if (tmp[i] == 'P')
					n = 97;
				if (tmp[i] == 'V')
					n = 99;
				if (tmp[i] == 'T')
					n = 101;
				if (tmp[i] == 'C')
					n = 103;
				if (tmp[i] == 'I' || tmp[i] == 'L')
					n = 113;
				if (tmp[i] == 'N')
					n = 114;
				if (tmp[i] == 'D')
					n = 115;
				if (tmp[i] == 'K' || tmp[i] == 'Q')
					n = 128;
				if (tmp[i] == 'E')
					n = 129;
				if (tmp[i] == 'M')
					n = 131;
				if (tmp[i] == 'H')
					n = 137;
				if (tmp[i] == 'F')
					n = 147;
				if (tmp[i] == 'R')
					n = 156;
				if (tmp[i] == 'Y')
					n = 163;
				if (tmp[i] == 'W')
					n = 186;
				sum += n;
			}
			arr[k] = sum;
			k++;
			sum = 0;
			tmp = "";
			if (j == pepide.length())
				break;
		}
	}
	arr = sort(arr);
	return arr;
}

string expand(string peptides, int k) {
	string amino = "GASPVTCINDKEMHFRYW";
	string str;

	if (peptides == "") 
		str = amino;

	for (int i = 0; i < peptides.length(); i = k + i) {
		string curr = peptides.substr(i,k);
		for (int j = 0; j < amino.length(); j++) {
			str += curr;
			str += amino.substr(j, 1);
		}
	}

	return str;
}

int parent_mass(vector<int> spectrum) {
	int mass = spectrum[spectrum.size() - 1];
	return mass;
}

int mass(string tmp) {
	int mass_tmp = 0;
	int n = 0;

	for (int i = 0; i < tmp.length(); i++) {
		if (tmp[i] == 'G')
			n = 57;
		if (tmp[i] == 'A')
			n = 71;
		if (tmp[i] == 'S')
			n = 87;
		if (tmp[i] == 'P')
			n = 97;
		if (tmp[i] == 'V')
			n = 99;
		if (tmp[i] == 'T')
			n = 101;
		if (tmp[i] == 'C')
			n = 103;
		if (tmp[i] == 'I' || tmp[i] == 'L')
			n = 113;
		if (tmp[i] == 'N')
			n = 114;
		if (tmp[i] == 'D')
			n = 115;
		if (tmp[i] == 'K' || tmp[i] == 'Q')
			n = 128;
		if (tmp[i] == 'E')
			n = 129;
		if (tmp[i] == 'M')
			n = 131;
		if (tmp[i] == 'H')
			n = 137;
		if (tmp[i] == 'F')
			n = 147;
		if (tmp[i] == 'R')
			n = 156;
		if (tmp[i] == 'Y')
			n = 163;
		if (tmp[i] == 'W')
			n = 186;
		mass_tmp += n;
	}
	return mass_tmp;
}

int is_consistent(string pepide, vector<int> spectrum) {

	int count = 0;
	int res = 0;
	int sum = 0;
	int n = 0;
	int size = pepide.length();
	vector<int> mas = { 0 };
	string tmp;

	for (int j = 1; j < pepide.length() + 1; j++) {
		for (int m = 0; m < pepide.length(); m++) {

			if (j == pepide.length())
				tmp = pepide;

			else if (m + j > pepide.length()) {
				break;
			}
			else tmp = pepide.substr(m, j);

			for (int i = 0; i < tmp.length(); i++) {
				if (tmp[i] == 'G')
					n = 57;
				if (tmp[i] == 'A')
					n = 71;
				if (tmp[i] == 'S')
					n = 87;
				if (tmp[i] == 'P')
					n = 97;
				if (tmp[i] == 'V')
					n = 99;
				if (tmp[i] == 'T')
					n = 101;
				if (tmp[i] == 'C')
					n = 103;
				if (tmp[i] == 'I' || tmp[i] == 'L')
					n = 113;
				if (tmp[i] == 'N')
					n = 114;
				if (tmp[i] == 'D')
					n = 115;
				if (tmp[i] == 'K' || tmp[i] == 'Q')
					n = 128;
				if (tmp[i] == 'E')
					n = 129;
				if (tmp[i] == 'M')
					n = 131;
				if (tmp[i] == 'H')
					n = 137;
				if (tmp[i] == 'F')
					n = 147;
				if (tmp[i] == 'R')
					n = 156;
				if (tmp[i] == 'Y')
					n = 163;
				if (tmp[i] == 'W')
					n = 186;
				sum += n;
			}
			mas.push_back(sum);
			sum = 0;
			tmp = "";
			if (j == pepide.length())
				break;
		}
	}
	mas = sort(mas);

	int z = 0;
	for (int i = 0; i < mas.size(); i++) {
		while (z < spectrum.size()) {
			if (mas[i] == spectrum[z]) {
				z++;
				count++;
				break;
			}
			else z++;
		}

	}
	if (count == mas.size())
		res = 1;
	else res = 0;
	return res;
}

int cyclospectrum(string pepide, vector<int> spectrum) {
	int count = 0;
	int res = 0;
	vector<int> mas = generating_theoretical_spectrum(pepide);

	for (int i = 0; i < mas.size(); i++) {
		if (mas[i] == spectrum[i])
			count++;
	}
	if (count == mas.size() && count == spectrum.size())
		res = 1;
	else res = 0;
	return res;
}

void cyclopeptide_sequencing(vector<int> spectrum) {
	string peptid = "";
	string peptides = "0";
	string tmp, text;
	int k = -1;
	int m = 0;
	while (!peptides.empty()) {
		if (k == -1)
			peptides = "";
		k++;
		m++;
		peptides = expand(peptides, k);
		int size_peptides = peptides.length();
		int curr = 0;
		for (int i = 0; i < size_peptides; i = i + m) {
			peptid = peptides.substr(curr,m);
			if (mass(peptid) == parent_mass(spectrum)) {
				if (cyclospectrum(peptid, spectrum) == 1) {

					tmp = "";
					tmp = peptid;

					int n = 0;
					int size = peptid.length()*(peptid.length() - 1) + 2;
					int *arr = new int[size];

					for (int i = 0; i < tmp.length(); i++) {
						if (tmp[i] == 'G')
							n = 57;
						if (tmp[i] == 'A')
							n = 71;
						if (tmp[i] == 'S')
							n = 87;
						if (tmp[i] == 'P')
							n = 97;
						if (tmp[i] == 'V')
							n = 99;
						if (tmp[i] == 'T')
							n = 101;
						if (tmp[i] == 'C')
							n = 103;
						if (tmp[i] == 'I' || tmp[i] == 'L')
							n = 113;
						if (tmp[i] == 'N')
							n = 114;
						if (tmp[i] == 'D')
							n = 115;
						if (tmp[i] == 'K' || tmp[i] == 'Q')
							n = 128;
						if (tmp[i] == 'E')
							n = 129;
						if (tmp[i] == 'M')
							n = 131;
						if (tmp[i] == 'H')
							n = 137;
						if (tmp[i] == 'F')
							n = 147;
						if (tmp[i] == 'R')
							n = 156;
						if (tmp[i] == 'Y')
							n = 163;
						if (tmp[i] == 'W')
							n = 186;
						arr[i] = n;
						if (i == tmp.length() - 1)
							cout << arr[i];
						else cout << arr[i] << "-";
					}
					cout << ' ';
				}
				peptides.erase(curr, m);
			}
			else if (is_consistent(peptid, spectrum) == 0) {
				peptides.erase(curr, m);
			}
			else curr = curr + m;
		}
	}
}

int cyclopeptide_scoring(string pepide, vector<int> spectrum) {

	vector<int> mas = generating_theoretical_spectrum(pepide);
	int count = 0;
	for (int i = 0; i < mas.size(); i++) {
		for (int z = 0; z < spectrum.size(); z++) {
			if (mas[i] == spectrum[z]) {
				spectrum[z] = 0;
				count++;
				break;
			}
		}
	}
	return count;
}


int main() {

	/*vector<int> all;
	int i = 0;
	int k = 0;
	while (i < 8) {
		cin >> k;

		if (cin.fail())
			break;

		all.push_back(k);
		i++;
	}
	cyclopeptide_sequencing(all);

	string p;
	cin >> p;
	vector<int> all;
	int i = 0;
	int k = 0;
	while (i < 1000) {
		cin >> k;

		if (cin.fail())
			break;

		all.push_back(k);
		i++;
	}
	int m = cyclopeptide_scoring(p, all);
	cout << m;*/

	return 0;
}