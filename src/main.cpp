#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>

using namespace std;

const int OPERATOR = 0;
const int CHAR = 1;

struct NFA {
    int s, t;
    vector<vector<vector<int>>> map;

    NFA(NFA &&n) noexcept : s(n.s), t(n.t), map(move(n.map)){
        n.s = n.t = 0;
    }
    NFA(int s, int t, vector<vector<vector<int>>> &&map) : s(s), t(t), map(map){}
};

struct DFA {
    vector<vector<int>> map;

    DFA(DFA &&n) noexcept : map(move(n.map)){}
    explicit DFA(vector<vector<int>> &&map) : map(map){}
};

struct FinalDFA {
    vector<vector<int>> map, end;

    FinalDFA(FinalDFA &&x) noexcept : map(move(x.map)), end(move(x.end)) {}
    FinalDFA(vector<vector<int>> &&map, vector<vector<int>> &&end) : map(map), end(end) {}
};

/*-----------函数声明-------------*/
int createCPP(ofstream&, const FinalDFA&, const vector<string>&);
int judge(char);
vector<pair<int, char>> preTreatment(const string&);
int optPriority(char);
vector<pair<int, char>> toPostfix(const string&);
NFA toNFA(const vector<pair<int, char>>&);
vector<int> eClosure(const vector<int>&, const vector<vector<vector<int>>>&);
DFA toDFA(const NFA&);
DFA optimizeDFA(const DFA&);
int DfsUpdate(vector<vector<int>>&, vector<int>&, int, int, const vector<vector<int>>&);
FinalDFA mergeDFAs(const vector<DFA>&);
/*-------------------------------*/


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Source file(*.l) is needed." << endl;
        return 1;
    }
    if (argc > 2)
        cout << "The redundant argument(s) will be ignored." << endl;
    string lexFile(argv[1]);
    ifstream ifs(lexFile, ifstream::in);
    if (!ifs) {
        cout << "No such file." << endl;
        return 1;
    }
    auto pos = lexFile.find_last_of('/') == -1 ? 0 : lexFile.find_last_of('/') + 1;
    auto n = lexFile.find_last_of('.') == -1 ? (lexFile.length() - pos - 1) : (lexFile.find_last_of('.') - pos);
    string name(lexFile, pos, n);
    ofstream ofs(name + ".yy.cpp", ofstream::out);

    vector<DFA> dfas;
    vector<string> tokenNames;
    while (ifs.good()) {
        string tokenName, tokenRE;
        ifs >> tokenName >> tokenRE;
        tokenNames.push_back(move(tokenName));
        dfas.push_back(optimizeDFA(toDFA(toNFA(toPostfix(tokenRE)))));
    }
    FinalDFA finalDFA = mergeDFAs(dfas);

    createCPP(ofs, finalDFA, tokenNames);

    ifs.close();
    ofs.close();
    return 0;
}

/**
 * 初始化输出文件
 * @param ofs
 * @return
 */
inline int createCPP(ofstream &ofs, const FinalDFA &dfa, const vector<string> &tokenNames) {
    ofs << "#include <iostream>\n"
            "#include <fstream>\n"
            "#include <vector>\n"
            "using namespace std;\n"
            "const vector<vector<int>> dfa{";
    auto &dMap = dfa.map;
    for (int i = 0; i < dMap.size(); ++i) {
        ofs << "{";
        for (int j = 0; j < dMap[i].size(); ++j) {
            ofs << dMap[i][j];
            if (j != dMap[i].size() - 1)
                ofs << ',';
        }
        ofs << '}';
        if (i != dMap.size() - 1)
            ofs << ',';
    }
    ofs << "};\nconst vector<vector<int>> dEnd{";
    auto &dEnd = dfa.end;
    for (int i = 0; i < dEnd.size(); ++i) {
        ofs << "{";
        for (int j = 0; j < dEnd[i].size(); ++j) {
            ofs << dEnd[i][j];
            if (j != dEnd[i].size() - 1)
                ofs << ',';
        }
        ofs << '}';
        if (i != dEnd.size() - 1)
            ofs << ',';
    }
    ofs << "};\nconst vector<string> tokenNames{";
    for (int i = 0; i < tokenNames.size(); ++i) {
        ofs << '\"' << tokenNames[i] << '\"';
        if (i != tokenNames.size() - 1)
            ofs << ',';
    }
    ofs << "};\nint main(int argc, char *argv[]) {\n"
            "    if (argc < 2) {\n"
            "        cout << \"Source file is needed.\" << endl;\n"
            "        return 1;\n"
            "    }\n"
            "    if (argc > 2)\n"
            "        cout << \"The redundant argument(s) will be ignored.\" << endl;\n"
            "    string src(argv[1]);\n"
            "    ifstream ifs(src, ifstream::in);\n"
            "    if (!ifs) {\n"
            "        cout << \"No such file.\" << endl;\n"
            "        return 1;\n"
            "    }\n"
            "    int state = 1;\n"
            "    char c = 0;\n"
            "    string word;\n"
            "    while (ifs.good() || c) {\n"
            "        if (!c) {\n"
            "            c = ifs.get();\n"
            "            if (c == EOF)\n"
            "                break;\n"
            "        }\n"
            "        if (!dfa[state][c]) {\n"
            "            if (word.empty() || dEnd[state].empty()) {\n"
            "                cout << \"Error.\" << endl;\n"
            "                return 1;\n"
            "            } else {\n"
            "                cout << '<' << tokenNames[dEnd[state][0]] << \", \\\"\" << word << \"\\\">\" << endl;\n"
            "                word = \"\";\n"
            "                state = 1;\n"
            "            }\n"
            "        } else {\n"
            "            string added;\n"
            "            switch (c) {\n"
            "                case '\\f':\n"
            "                    added = \"\\\\f\";\n"
            "                    break;\n"
            "                case '\\n':\n"
            "                    added = \"\\\\n\";\n"
            "                    break;\n"
            "                case '\\r':\n"
            "                    added = \"\\\\r\";\n"
            "                    break;\n"
            "                case '\\t':\n"
            "                    added = \"\\\\t\";\n"
            "                    break;\n"
            "                case '\\v':\n"
            "                    added = \"\\\\v\";\n"
            "                    break;\n"
            "                default:\n"
            "                    added = c;\n"
            "                    break;\n"
            "            }\n"
            "            word += added;\n"
            "            state = dfa[state][c];\n"
            "            c = 0;\n"
            "        }\n"
            "    }\n"
            "    if (!word.empty() && !dEnd[state].empty())\n"
            "        cout << '<' << tokenNames[dEnd[state][0]] << \", \\\"\" << word << \"\\\">\" << endl;\n"
            "    else {\n"
            "        cout << \"Error.\" << endl;\n"
            "        return 1;\n"
            "    }\n"
            "    return 0;\n"
            "}";
    return 0;
}

/**
 * 判断字符类型
 * @param c
 * @return
 */
inline int judge(char c) {
    switch (c) {
        case '(':
        case ')':
        case '[':
        case ']':
        case '*':
        case '|':
        case '-':
            return OPERATOR;
        default:
            return CHAR;
    }
}

/**
 * 规范化正规表达式
 * @param re0
 * @return
 */
inline vector<pair<int, char>> preTreatment(const string &re0) {
    // 替换转义字符
    vector<pair<int, char>> re1;
    for (int i = 0; i < re0.length(); ++i) {
        if (re0[i] == '\\') {
            switch (re0[i + 1]) {
                case 'd':
                    re1.emplace_back(make_pair(OPERATOR, '['));
                    re1.emplace_back(make_pair(CHAR, '0'));
                    re1.emplace_back(make_pair(OPERATOR, '-'));
                    re1.emplace_back(make_pair(CHAR, '9'));
                    re1.emplace_back(make_pair(OPERATOR, ']'));
                    break;
                case 's':
                    re1.emplace_back(make_pair(OPERATOR, '['));
                    re1.emplace_back(make_pair(CHAR, ' '));
                    re1.emplace_back(make_pair(CHAR, '\f'));
                    re1.emplace_back(make_pair(CHAR, '\n'));
                    re1.emplace_back(make_pair(CHAR, '\r'));
                    re1.emplace_back(make_pair(CHAR, '\t'));
                    re1.emplace_back(make_pair(CHAR, '\v'));
                    re1.emplace_back(make_pair(OPERATOR, ']'));
                    break;
                default:
                    re1.emplace_back(CHAR, re0[i + 1]);
                    break;
            }
            ++i;
        } else {
            re1.emplace_back(judge(re0[i]), re0[i]);
        }
    }
    // 增加连接符，替换[]
    vector<pair<int, char>> newRE;
    stack<char> bracket;
    if (re1[0].first == OPERATOR && (re1[0].second == '(' || re1[0].second == '[')) {
        bracket.push(re1[0].second);
        newRE.emplace_back(OPERATOR, '(');
    } else
        newRE.push_back(re1[0]);
    for (int i = 1, j = 0; i < re1.size(); ++i) {
        if ((re1[i].first == CHAR || re1[i].second == '(' || re1[i].second == '[')
                && (newRE[j].first == CHAR || newRE[j].second == ')' || newRE[j].second == '*')) {
            newRE.emplace_back(OPERATOR, (bracket.empty() || bracket.top() == '(') ? '.' : '|');
            ++j;
        }
        if (re1[i].first == OPERATOR) {
            if (re1[i].second == '[' || re1[i].second == '(')
                bracket.push(re1[i].second);
            if (re1[i].second == ']' || re1[i].second == ')')
                bracket.pop();
            if (re1[i].second == '[') {
                newRE.emplace_back(OPERATOR, '(');
                ++j;
                continue;
            }
            if (re1[i].second == ']') {
                newRE.emplace_back(OPERATOR, ')');
                ++j;
                continue;
            }
            if (re1[i].second == '-') {
                newRE[j].second = '(';
                newRE[j].first = OPERATOR;
                newRE.emplace_back(CHAR, re1[i - 1].second);
                ++j;
                for (char t = (char)(re1[i - 1].second + 1); t <= re1[i + 1].second; ++t) {
                    newRE.emplace_back(OPERATOR, '|');
                    newRE.emplace_back(CHAR, t);
                    j += 2;
                }
                newRE.emplace_back(OPERATOR, ')');
                ++j;
                ++i;
            } else {
                newRE.emplace_back(OPERATOR, re1[i].second);
                ++j;
            }
        } else {
            newRE.emplace_back(CHAR, re1[i].second);
            ++j;
        }
    }
    return newRE;
}

/**
 * 返回运算符的优先级
 * @param c
 * @return
 */
inline int optPriority(char c) {
    switch (c) {
        case '|':
            return 1;
        case '.':
            return 2;
        case '*':
            return 3;
        case '(':
        case ')':
            return 4;
        default:
            return 0;
    }
}

/**
 * 将正则表达式的中缀形式转为后缀形式
 * @param origin
 * @return
 */
inline vector<pair<int, char>> toPostfix(const string &origin) {
    auto newInfix = preTreatment(origin);
    vector<pair<int, char>> postfix;
    stack<char> opt;
    for (auto c: newInfix) {
        if (c.first == CHAR)
            postfix.emplace_back(CHAR, c.second);
        else {
            switch (c.second) {
                case '(':
                    opt.push(c.second);
                    break;
                case ')':
                    while (!opt.empty() && opt.top() != '(') {
                        postfix.emplace_back(OPERATOR, opt.top());
                        opt.pop();
                    }
                    opt.pop();
                    break;
                default:
                    while (!opt.empty() && opt.top() != '(' && optPriority(opt.top()) >= optPriority(c.second)) {
                        postfix.emplace_back(OPERATOR, opt.top());
                        opt.pop();
                    }
                    opt.push(c.second);
            }
        }
    }
    while (!opt.empty()) {
        postfix.emplace_back(OPERATOR, opt.top());
        opt.pop();
    }
    return postfix;
}

/**
 * re --> NFA
 * @param re
 * @return
 */
inline NFA toNFA(const vector<pair<int, char>> &re) {
    vector<vector<vector<int>>> nfa(1, vector<vector<int>>(128));
    stack<pair<int, int>> subNFA;
    for (auto i: re) {
        if (i.first == CHAR) {
            int s = (int)nfa.size(), t = s + 1;
            nfa.emplace_back(vector<vector<int>>(128));
            nfa.emplace_back(vector<vector<int>>(128));
            nfa[s][i.second].emplace_back(t);
            subNFA.emplace(s, t);
        } else {
            if (i.second == '.') {
                pair<int, int> right = move(subNFA.top());
                subNFA.pop();
                pair<int, int> left = move(subNFA.top());
                subNFA.pop();
                nfa[left.second][0].emplace_back(right.first);
                subNFA.emplace(left.first, right.second);
            } else if (i.second == '|') {
                pair<int, int> e1 = move(subNFA.top());
                subNFA.pop();
                pair<int, int> e2 = move(subNFA.top());
                subNFA.pop();
                int newLeft = (int)(nfa.size()), newRight = newLeft + 1;
                nfa.emplace_back(vector<vector<int>>(128));
                nfa.emplace_back(vector<vector<int>>(128));
                nfa[newLeft][0].emplace_back(e1.first);
                nfa[newLeft][0].emplace_back(e2.first);
                nfa[e1.second][0].emplace_back(newRight);
                nfa[e2.second][0].emplace_back(newRight);
                subNFA.emplace(newLeft, newRight);
            } else if (i.second == '*') {
                pair<int, int> cir = move(subNFA.top());
                subNFA.pop();
                int newState = (int)nfa.size();
                nfa.emplace_back(vector<vector<int>>(128));
                nfa[newState][0].emplace_back(cir.first);
                nfa[cir.second][0].emplace_back(newState);
                subNFA.emplace(newState, newState);
            }
        }
    }
    return NFA(subNFA.top().first, subNFA.top().second, move(nfa));
}

/**
 * ε-闭包
 * @param states
 * @param map
 * @return
 */
inline vector<int> eClosure(const vector<int> &states, const vector<vector<vector<int>>> &map) {
    vector<int> result;
    queue<int> q;
    vector<bool> flag(map.size());
    for (int node: states) {
        q.push(node);
        flag[node] = true;
    }
    while (!q.empty()) {
        int front = q.front();
        q.pop();
        result.push_back(front);
        for (int to: map[front][0])
            if (!flag[to])
                q.push(to);
    }
    return result;
}

/**
 * NFA --> DFA
 * @param nfa
 * @return
 */
inline DFA toDFA(const NFA &nfa) {
    // dfaMap[i][0] 表示新状态I(i)是否为终态
    vector<vector<int>> dfaMap(1, vector<int>());
    vector<vector<int>> newStates(1, vector<int>());
    vector<int> firstState = eClosure(vector<int>{nfa.s}, nfa.map);
    sort(firstState.begin(), firstState.end());
    newStates.push_back(firstState);
    dfaMap.emplace_back(vector<int>(128));
    for (int state: firstState)
        if (state == nfa.t) {
            dfaMap[1][0] = 1;
            break;
        }
    for (int i = 1; i < newStates.size(); ++i) {
        for (int j = 1; j < 128; ++j) {
            vector<int> temp;
            for (int k = 0; k < newStates[i].size(); ++k) {
                if (!nfa.map[newStates[i][k]][j].empty()) {
                    bool exist = false;
                    for (int t: temp)
                        if (t == nfa.map[newStates[i][k]][j][0]) {
                            exist = true;
                            break;
                        }
                    if (!exist)
                        temp.push_back(nfa.map[newStates[i][k]][j][0]);
                }
            }
            if (temp.empty())
                continue;
            temp = eClosure(temp, nfa.map);
            sort(temp.begin(), temp.end());
            int index = 0;
            for (int ii = 1; ii < newStates.size(); ++ii) {
                if (newStates[ii].size() != temp.size())
                    continue;
                bool same = true;
                for (int jj = 0; jj < temp.size(); ++jj)
                    if (newStates[ii][jj] != temp[jj]) {
                        same = false;
                        break;
                    }
                if (same) {
                    index = ii;
                    break;
                }
            }
            if (index)
                dfaMap[i][j] = index;
            else {
                dfaMap[i][j] = (int)(newStates.size());
                int end = 0;
                for (int state: temp)
                    if (state == nfa.t) {
                        end = 1;
                        break;
                    }
                newStates.emplace_back(temp);
                dfaMap.emplace_back(vector<int>(128));
                dfaMap[newStates.size() - 1][0] = end;
            }
        }
    }
    return DFA(move(dfaMap));
}

/**
 * 优化DFA，减少状态数
 * @param dfa
 * @return
 */
inline DFA optimizeDFA(const DFA &dfa) {
    auto &dMap = dfa.map;
    vector<vector<int>> newStates(3);
    vector<int> mark(dfa.map.size());
    for (int i = 1; i < dMap.size(); ++i)
        if (!dMap[i][0]) {
            newStates[1].push_back(i);
            mark[i] = 1;
        }
        else {
            newStates[2].push_back(i);
            mark[i] = 2;
        }
    mark[0] = 3;
    bool updated = true;
    while (updated) {
        updated = false;
        int oldSize = (int)newStates.size();
        for (int i = 1; i < oldSize; ++i)
            updated = updated || DfsUpdate(newStates, mark, i, 1, dfa.map);
    }
    vector<vector<int>> newMap(newStates.size(), vector<int>(128));
    mark[0] = 0;
    for (int i = 1; i < newStates.size(); ++i) {
        for (int state: newStates[i]) {
            for (int c = 1; c < 128; ++c)
                newMap[i][c] = mark[dMap[state][c]];
            newMap[i][0] = newMap[i][0] || dMap[state][0];
        }
    }
    return DFA(move(newMap));
}

/**
 * 深度优先，划分不等价的状态集
 * @param newStates
 * @param mark
 * @param index
 * @param from
 * @param dMap
 * @return
 */
int DfsUpdate(vector<vector<int>> &newStates, vector<int> &mark, int index, int transferChar, const vector<vector<int>> &dMap) {
    if (transferChar > 127 || newStates[index].size() == 1)
        return 0;
    int oldSize = (int)newStates.size();
    vector<int> &updating = newStates[index];
    vector<pair<int, vector<int>>> temp;
    for (int state: updating) {
        int to = dMap[state][transferChar];
        bool exist = false;
        for (auto &_i: temp)
            if (mark[to] == _i.first) {
                _i.second.push_back(state);
                exist = true;
                break;
            }
        if (!exist)
            temp.emplace_back(mark[to], vector<int>{state});
    }
    if (temp.size() > 1) {
        newStates[index] = vector<int>();
        for (int ii: temp[0].second)
            newStates[index].push_back(ii);
        for (int ii = 1; ii < temp.size(); ++ii) {
            newStates.emplace_back(vector<int>());
            for (int _state: temp[ii].second) {
                newStates[newStates.size() - 1].push_back(_state);
                mark[_state] = mark[0];
            }
            ++mark[0];
        }
        int newSize = (int)newStates.size();
        DfsUpdate(newStates, mark, index, transferChar + 1, dMap);
        for (int i = oldSize; i < newSize; ++i)
            DfsUpdate(newStates, mark, i, transferChar + 1, dMap);
        return 1;
    }
    return DfsUpdate(newStates, mark, index, transferChar + 1, dMap);
}

/**
 * 合并多个DFA
 * @return
 */
inline FinalDFA mergeDFAs(const vector<DFA> &dfas) {
    vector<vector<pair<int, int>>> newStates(2);
    vector<vector<int>> newMap(2, vector<int>(128));
    vector<vector<int>> end(2);
    for (int i = 0; i < dfas.size(); ++i)
        newStates[1].emplace_back(i, 1);
    for (int i = 1; i < newStates.size(); ++i)
        for (int c = 1; c < 128; ++c) {
            vector<pair<int, int>> temp;
            for (auto state: newStates[i])
                if (dfas[state.first].map[state.second][c])
                    temp.emplace_back(state.first, dfas[state.first].map[state.second][c]);
            if (temp.empty())
                continue;
            sort(temp.begin(), temp.end());
            int index = 0;
            for (int ii = 1; ii < newStates.size(); ++ii) {
                if (newStates[ii].size() != temp.size())
                    continue;
                bool same = true;
                for (int jj = 0; jj < temp.size(); ++jj)
                    if (newStates[ii][jj] != temp[jj]) {
                        same = false;
                        break;
                    }
                if (same) {
                    index = ii;
                    break;
                }
            }
            if (index)
                newMap[i][c] = index;
            else {
                newMap.emplace_back(vector<int>(128));
                newMap[i][c] = (int)newStates.size();
                end.emplace_back(vector<int>());
                for (auto &ss: temp)
                    if (dfas[ss.first].map[ss.second][0])
                        end[newStates.size()].emplace_back(ss.first);
                newStates.emplace_back(temp);
            }
        }
    return FinalDFA(move(newMap), move(end));
}