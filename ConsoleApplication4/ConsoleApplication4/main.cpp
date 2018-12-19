#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
using namespace std;

vector<string> permute(vector<string> tmp) {
	string s = "ACGT";
	vector<string> res, temp;
	temp = tmp;
	for (int i = 0; i < s.length(); i++)
	{
		for (int j = 0; j < tmp.size(); j++)
		{
			res.push_back(tmp[j] + s[i]);
		}
		tmp = temp;
	}
	return res;
}
vector<string> f(string tmp,vector<string> tmp1, int d) {
	vector<string> res1;
	res1.push_back(tmp);
	for (int m = 0; m < tmp1.size(); m++) {
		string v4 = tmp1[m];
		int count = 0;
		for (int mm = 0; mm < v4.length(); mm++)
			if (v4[mm] != tmp[mm])
				count++;
		if (count == d)
			res1.push_back(v4);
	}
	return res1;
}
void motif_enumeration(vector<string> dna, int k, int d, int n) {
	vector<string> v1, v2, res;
	vector<int> res_count;
	string s = dna[0];
	string s1;
	int count;
	set<string> mS;
	string tmp, tmp1, tmp2, tmp2_k;
	v1.push_back("A");
	v1.push_back("C");
	v1.push_back("G");
	v1.push_back("T");

	int j = 1;
	while (j < k)
	{
		v2 = permute(v1);
		v1 = v2;
		j++;
	}
	vector<string> res1;
	vector<string> res2;
	string motif;
	for (int z = 0; z < s.length() - k + 1; z++) {
		tmp = s.substr(z, k);
		res1.clear();
		res1 = f(tmp, v2, d);
		for (int i = 0; i < res1.size(); i++) {
			motif = "";
			motif = res1[i];
			count = 0;
			for (int num = 1; num < n; num++) {
				tmp2 = dna[num];
				for (int mm = 0; mm < tmp2.length()- k + 1; mm++) {
					tmp2_k = tmp2.substr(mm, k);

					int chislo_mut = 0;
					for (int m = 0; m < motif.length(); m++) {
						if (motif[m] != tmp2_k[m])
							chislo_mut++;
					}
					if (d == chislo_mut || chislo_mut == 0) {
						count++;
						break;
					}
				}
			}
			if (count >= n - 1)
				mS.insert(motif);
		}
	
	}
	
	for (std::set<string>::iterator it = mS.begin(); it != mS.end(); ++it)
	{
		std::cout << *it << ' ';
	}
}

int d(vector<string> dna, string motifs) {
	int size = dna.size();
	string tmp;
	int sum = 0;
	int k = 0;
	int min, curr;
	int m = 0;
	int n = 0;

	for (int i = 0; i < size; i++) {
		tmp = dna[i];
		if (motifs.length() == tmp.length())
			m = 1;
		for (int j = 0; j < tmp.length(); j++) {
			if (m == 1)
				break;
			k = 0;
			curr = j;

			while (n < motifs.length() && j < tmp.length()) {
				if (motifs[n] != tmp[j])
					k = k + 1;
				j++;
				n++;
			}
			m = m - 1;
			j = curr;

			if (n != motifs.length()) {
				n = 0;
				break;
			}
			n = 0;
			if (j == 0)
				min = k;
			else {
				if (min > k)
					min = k;
			}
		}
		sum += min;
	}
	return sum;
}
string median_string(vector<string> dna, int k){

	string res = "";
	vector<string> motifs;

	motifs.push_back("A");
	motifs.push_back("C");
	motifs.push_back("G");
	motifs.push_back("T");

	for (int i = 0; i < k -1; i++)
		motifs = permute(motifs);

	int size = motifs.size();
	int dist = dna[0].length();

	for (int i = 0; i < size; i++)
	{
		int tmp = d(dna, motifs[i]);
		if (tmp < dist)
		{
			res = motifs[i];
			dist = tmp;
		}
	}
	return res;
}

int main() {
	/*int k, d;
	vector<string> dna;
	string s;
	cin >> k;
	cin >> d;

	int n = 0;
	while (!cin.eof())
	{
		cin >> s;
		n++;
		dna.push_back(s);
	}
	motif_enumeration(dna,k,d,n);*/

	int m;
	vector<string> dna1;
	string s1;
	while (!cin.eof())
	{
		cin >> s1;
		dna1.push_back(s1);
	}
	string s2;
	s2 = median_string(dna1, m);
	for (int i = 0; i < s2.length(); i++)
		cout << s2[i];
	cout << endl;
	return 0;
}