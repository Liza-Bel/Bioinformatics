#include <iostream>
#include <string>
#include <map>
#include <algorithm>
using namespace std;

size_t pattern_count(string pattern, string genome ) {
	size_t count = 0;
	size_t pos = pattern.find(genome);

	while (pos != std::string::npos) {
		++count;
		++pos;
		pos = pattern.find(genome, pos);
	}

	return count;
}

void frequent_words(string text, int k) { 
	map<string, int> mymap;
	int max = 0;

	for (int i = 0; i < text.length(); i++) {
		string tmp = text.substr(i, k);
		mymap[tmp]++;
		if (max < mymap[tmp])
			max = mymap[tmp];
	}

	for (auto i = mymap.begin(); i != mymap.end(); i++) {
		if (i->second == max)
			cout << i->first << endl;
	}
}

string revers_complement(string pattern) {
	reverse(pattern.begin(), pattern.end());
	for (int i = 0; i < pattern.length() ; i++) {
		if (pattern[i] == 'A') {
			pattern[i] = 'T';
		}
		else if (pattern[i] == 'T') {
			pattern[i] = 'A';
		}
		else if (pattern[i] == 'C') {
			pattern[i] = 'G';
		}
		else  {
			pattern[i] = 'C';
		}
	}
	return pattern;
}


int main() {
	cout << "Pattern Count" << endl;
	string a, b;
	cout << "Pattern: ";
	cin >> a;
	cout << "Genome: ";
	cin >> b;
	size_t count = pattern_count(a, b);
	cout << "Count: " << count << endl;
	cout << endl;

	cout << "Frequent Words" << endl;
	string c;
	int k;
	cout << "Pattern: ";
	cin >> c;
	cout << "K: ";
	cin >> k;
	frequent_words(c, k);
	cout << "Res: " << endl;
	cout << endl;

	cout << "Revers Complement" << endl;
	cout << "Pattern:";
	string d;
	cin >> d;
	string res = revers_complement(d);
	cout << "Res: ";
	cout << res;

	return 0;
}