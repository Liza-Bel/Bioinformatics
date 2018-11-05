#include <iostream>
#include <string>
#include <map>
#include <algorithm>
using namespace std;

string protein_translatoin(string pattern) {
	string tmp;
	int j = 0;
	string n;
	for (int i = 0; i < pattern.length(); i = i + 3) {
		string tmp = pattern.substr(i, 3);

		if (tmp == "UGA" || tmp == "UAG" || tmp == "UAA")
			continue;

		if (tmp == "AAU" || tmp == "AAC")
			pattern[j] = 'N';
		else if (tmp == "AAA" || tmp == "AAG")
			pattern[j] = 'K';
		else if (tmp == "ACU" || tmp == "ACC" || tmp == "ACA" || tmp == "ACG")
			pattern[j] = 'T';
		else if (tmp == "AGU" || tmp == "AGC")
			pattern[j] = 'S';
		else if (tmp == "AGA" || tmp == "AGG")
			pattern[j] = 'R';
		else if (tmp == "AUU" || tmp == "AUC" || tmp == "AUA")
			pattern[j] = 'I';
		else if (tmp == "AUG")
			pattern[j] = 'M';

		else if (tmp == "CAU" || tmp == "CAC")
			pattern[j] = 'H';
		else if (tmp == "CAA" || tmp == "CAG")
			pattern[j] = 'Q';
		else if (tmp == "CCU" || tmp == "CCC" || tmp == "CCA" || tmp == "CCG")
			pattern[j] = 'P';
		else if (tmp == "CGU" || tmp == "CGC" || tmp == "CGA" || tmp == "CGG")
			pattern[j] = 'R';
		else if (tmp == "CUU" || tmp == "CUC" || tmp == "CUA" || tmp == "CUG")
			pattern[j] = 'L';

		else if (tmp == "GAU" || tmp == "GAC")
			pattern[j] = 'D';
		else if (tmp == "GAA" || tmp == "GAG")
			pattern[j] = 'E';
		else if (tmp == "GCU" || tmp == "GCC" || tmp == "GCA" || tmp == "GCG")
			pattern[j] = 'A';
		else if (tmp == "GGU" || tmp == "GGC" || tmp == "GGA" || tmp == "GGG")
			pattern[j] = 'G';
		else if (tmp == "GUU" || tmp == "GUC" || tmp == "GUA" || tmp == "GUG")
			pattern[j] = 'V';

		else if (tmp == "UUG" || tmp == "UUA")
			pattern[j] = 'L';
		else if (tmp == "UUC" || tmp == "UUU")
			pattern[j] = 'F';
		else if (tmp == "UGG")
			pattern[j] = 'W';
		else if (tmp == "UGC" || tmp == "UGU")
			pattern[j] = 'C';
		else if (tmp == "UCG" || tmp == "UCA" || tmp == "UCC" || tmp == "UCU")
			pattern[j] = 'S';
		else if (tmp == "UAC" || tmp == "UAU")
			pattern[j] = 'Y';

		j++;
		
	}

	return pattern = pattern.substr(0, j);
}

string revers_complement(string pattern) {
	reverse(pattern.begin(), pattern.end());
	for (int i = 0; i < pattern.length(); i++) {
		if (pattern[i] == 'A') {
			pattern[i] = 'U';
		}
		else if (pattern[i] == 'U') {
			pattern[i] = 'A';
		}
		else if (pattern[i] == 'C') {
			pattern[i] = 'G';
		}
		else {
			pattern[i] = 'C';
		}
	}
	return pattern;
}

void peptide_encoding(string pattern, string amino) {

	for (int i = 0; i < pattern.length(); i++) {
		if (pattern[i] == 'T')
			pattern[i] = 'U';
	}

	string tmp, s, str;

	for (int i = 0; i < pattern.length(); i++) {
		tmp = pattern.substr(i, amino.length() * 3);

		s = protein_translatoin(tmp);
		str = revers_complement(tmp);
		str = protein_translatoin(str);

		if (amino == s || amino == str) {

			for (int i = 0; i < tmp.length(); i++) {
				if (tmp[i] == 'U')
					tmp[i] = 'T';
			}
			cout << tmp << endl;
		}
	}
}

long long int subpeptides_count_problem(long long int n) {
	return n*(n - 1);
}

void generating_theoretical_spectrum(string pepide) {
	int n = 0;
	int size = pepide.length()*(pepide.length() - 1) + 2;
	int *arr = new int[size];

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
				tmp = pepide.substr(m , pepide.length() - m);
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
	sort(arr,arr + size);
	for (int i = 0; i < size; i++) {
		cout << arr[i] << ' ';
	}
	delete arr;
}

long long int counting_peptides_with_given_mass(int n) {
	int arr[] = { 57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186 };
	map<int, long long int> mymap;
	mymap[0] = 1;

	for (int i = 57; i < n + 1; i++)
	{
		for (int j = 0; j < 18; j++)
			if (mymap[i - arr[j]] > 0)
				mymap[i] += mymap[i - arr[j]];
	}
	return mymap[n];
}

int main(){
	/*cout << "Protein Translatoin" << endl;
	string s,n;
	cin >> s;
	n = protein_translatoin(s);
	cout << n << endl;*/


	/*cout << "Peptide Encoding" << endl;
	string text, amino;
	cin >> text;
	cin >> amino;
	peptide_encoding(text, amino);*/
	
	/*cout << "Subpeptides Count Problem" << endl;
	long long int n1;
	cin >> n1;
	long long int k = subpeptides_count_problem(n1);
	cout << k << endl;

	cout << "Generating Theoretical Spectrum" << endl;
	string text1;
	cin >> text1;
	generating_theoretical_spectrum(text1);*/

	cout << "Counting Peptides With Given Mass" << endl;
	int m;
	cin >> m;
	long long int m1 = counting_peptides_with_given_mass(m);
	cout << m1 << endl;

	return 0;
}