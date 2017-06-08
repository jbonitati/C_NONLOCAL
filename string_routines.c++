void word_check (const class basic_string<char> words)
{
  class basic_string<char> word_to_explain;
  cin>>word_to_explain;

  if (word_to_explain == words) return;

  cout<<"Bad word : "<<word_to_explain<<endl;
  exit (1);
}

void word_check2 (const class basic_string<char> what_is_read,const class basic_string<char> words)
{

  if (what_is_read == words) return;

  cout<<"Bad word : "<<what_is_read<<endl;
  exit (1);
}

int determine_n(class basic_string<char> shell)
{
  class basic_string<char> n_string;

  for (int i = 0 ; shell[i] <= '9' ; i++)
    n_string += shell[i];

  return (atoi (shell.c_str ()));
}

int determine_l (class basic_string<char> shell)
{
  shell += '\n';

  int i = 0;
 
  while (true)
    switch (shell[i])
    {
    case 's': return 0;
    case 'p': return 1;
    case 'd': return 2;
    case 'f': return 3;
    case 'g': return 4;
    case 'h': return 5;
    case 'i': return 6;
    case 'j': return 7;
    case 'k': return 8;
    case 'l': return 9;
    case 'm': return 10;
    case 'n': return 11;
    case 'o': return 12;

    case '\n':
      cout<<"No orbital quantum number."<<endl;
      exit (1);

    default:
      i++;
    }
}



double determine_j (class basic_string<char> shell)
{
  shell += '\n';

  class basic_string<char> j_string;
  int i = 0;

  while (shell[i] <= '9') 
    i++;

  i++;

  while (shell[i] != '/')
  {
    j_string += shell[i];
    i++;

    if (shell[i] == '\n') return (atof (j_string.c_str ()));
  }

  return (atof (j_string.c_str ())/2.0);
}



int determine_J (class basic_string<char> eigenstate)
{
  class basic_string<char> J_string;

  for (int i = 0 ; ((eigenstate[i] != '+') && (eigenstate[i] != '-')) ; i++)
    J_string += eigenstate[i];

  return (atoi (J_string.c_str ()));
}


int determine_Parity (class basic_string<char> eigenstate)
{
  eigenstate += '\n';

  int i = 0;

  while (1)
    switch (eigenstate[i])
    {
    case '+':
      return 1;
      
    case '-':
      return -1;
    
    case '\n':
      cout<<"No parity."<<endl;
      exit (1);

    default:
      i++;
    }
}

