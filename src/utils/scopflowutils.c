#include <string.h>
#include <stdlib.h>
#include <stdio.h>
/**
 * Find first instance of character in a string after start
 * @param token character being searched for
 * @param string character string being searched
 * @param start index of first character being searched
 * @return index of first character corresponding to token after start
 * character. Return -1 if no character found
 */
int find_char(char token, char *string, int start)
{
  int len = strlen(string);
  int i;
  int ret = -1;
  for (i=start; i<len; i++) {
    if (string[i] == token) {
      ret = i;
      return ret;
    }
  }
  return ret;
}

/**
 * Find first instance of character in a string after start that does not match
 * token
 * @param token character being checked
 * @param string character string being searched
 * @param start index of first character being checked
 * @return index of first character that does not match token after start
 * character. Return -1 if no character does not match token
 */
int find_first_not_of(char token, char *string, int start)
{
  int len = strlen(string);
  int i;
  int ret = -1;
  for (i=start; i<len; i++) {
    if (string[i] != token) {
      ret = i;
      return ret;
    }
  }
  return ret;
}

/**
 * Clean up 2 character tags so that single quotes are removed and single
 * character tags are right-justified. These tags can be delimited by a
 * pair of single quotes, a pair of double quotes, or no quotes
 * @param string original string before reformatting and returns proper 2
 * character tag
 */
void clean2Char(char *string)
{
  int len = strlen(string);
  char *tag = (char*)malloc(sizeof(char)*len);;
  int ntok1, ntok2, sngl_qt, no_qt, i;
  /* Find and remove single or double quotes */
  ntok1 = find_char('\'',string,0);
  sngl_qt = 1;
  no_qt = 0;
  /* if no single quote found, then assume double quote or no quote */
  if (ntok1 == -1) {
    ntok1 = find_char('\"',string,0);
    /* if no double quote found then assume no quote */
    if (ntok1 == -1) {
      ntok1 = find_first_not_of(' ',string,0);
      no_qt = 1;
    } else {
      sngl_qt = 0;
    }
  }
  if (sngl_qt) {
    ntok1 = find_first_not_of('\'',string,ntok1);
    ntok2 = find_char('\'',string,ntok1);
  } else if (no_qt) {
    ntok2 = find_char(' ',string,ntok1);
  } else {
    ntok1 = find_first_not_of('\"',string,ntok1);
    ntok2 = find_char('\"',string,ntok1);
  }
  if (ntok2 == -1) ntok2 = len;
  for (i=0; i<ntok2-ntok1; i++) {
    string[i] = string[i+ntok1];
  }
  string[ntok2-ntok1] = '\0';
  /* get rid of white space */
  len = strlen(string);
  ntok1 = find_first_not_of(' ',string,0);
  ntok2 = find_char(' ',string,ntok1);
  if (ntok2 == -1) ntok2 = len;
  for (i=0; i<ntok2-ntok1; i++) {
    string[i] = string[i+ntok1];
  }
  string[ntok2-ntok1] = '\0';
  len = strlen(string);
  if (len == 1) {
    string[1] = ' ';
  }
  string[2] = '\0';
  return;
}

/**
 * Tokenize a string on blanks and return an array of strings
 * Blanks within single or double quotes are ignored. Note that array returned
 * by program must be free'd by calling program.
 * @param str string of tokens separated by blank characters
 * @param numtok number of tokens found
 * @param maxtokens maximumn number of tokens that are returned
 * @param maxchar maximumn number of characters in each token
 * @return array of tokens
 */
char** blankTokenizer(const char *str, int *numtok, int maxtokens, int maxchar)
{
  char **ret = (char**)malloc(sizeof(char*)*maxtokens);
  int i;
  int slen = strlen(str);
  char* strcpy = (char*)malloc(sizeof(char)*(slen+1));
  int ntok1, ntok2;

  for (i=0; i<maxtokens; i++) {
    ret[i] = (char*)malloc(sizeof(char)*maxchar);
  }
  /* ignore newline character if it is present */
  if (str[slen-1] == '\n') slen--;
  /* Replace any tabs with blank spaces */
  for (i=0; i<slen; i++) {
    if (strcpy[i] == '\t') {
      strcpy[i] = ' ';
    } else {
      strcpy[i] = str[i];
    }
  }
  ntok1 = find_first_not_of(' ',strcpy,0);
  if (strcpy[ntok1] == '\'') {
    ntok2 = find_char('\'',strcpy,ntok1+1);
    ntok2++;
  } else if (strcpy[ntok1] == '\"') {
    ntok2 = find_char('\"',strcpy,ntok1+1);
    ntok2++;
  } else if (ntok1 != -1) {
    ntok2 = find_char(' ',strcpy,ntok1);
    if (ntok2 == -1) ntok2 = slen;
  } else {
    *numtok = 0;
    return ret;
  }
  for (i=0; i<ntok2-ntok1; i++) {
    (ret[0])[i] = strcpy[i+ntok1];
  }
  (ret[0])[ntok2-ntok1] = '\0';
  *numtok = 1;
  while (ntok2 < slen-1 && ntok1 != -1) {
    ntok1 = find_first_not_of(' ',strcpy,ntok2);
    if (strcpy[ntok1] == '\'') {
      ntok2 = find_char('\'',strcpy,ntok1+1);
      if (ntok2 != -1) {
        ntok2++;
      } else {
        ntok2 = slen;
      }
    } else if (strcpy[ntok1] == '\"') {
      ntok2 = find_char('\"',strcpy,ntok1+1);
      if (ntok2 != -1) {
        ntok2++;
      } else {
        ntok2 = slen;
      }
    } else if (ntok1 != -1) {
      ntok2 = find_char(' ',strcpy,ntok1);
      if (ntok2 == -1) ntok2 = slen;
    } 
    if (ntok2 != -1) {
      for (i=0; i<ntok2-ntok1; i++) {
        (ret[*numtok])[i] = strcpy[i+ntok1];
      }
      (ret[*numtok])[ntok2-ntok1] = '\0';
      (*numtok)++;
    }
  }
  return ret;
}

/**
 * test program for tokenizer and clear2Char function
 *
int main(int argc, char **argv) {
  const char *string = "Flowers are natures \"best form\" of 3 \'beers\' return\n";
  int maxtokens = 10;
  int maxchar = 32;
  int numtok;
  char tag[32];
  char **tokens;
  int i;
  tokens = blankTokenizer(string, &numtok, maxtokens, maxchar);
  for (i=0; i<numtok; i++) {
    printf("token[%d]: (%s)\n",i,tokens[i]);
  }
  for (i=0; i<maxtokens; i++) {
    free(tokens[i]);
  }
  free(tokens);
  strncpy(tag,"\" C\"",4);
  clean2Char(tag);
  printf("clean2Char test1: (%s)\n",tag);
  strncpy(tag,"\' C\'",4);
  clean2Char(tag);
  printf("clean2Char test2: (%s)\n",tag);
  strncpy(tag," C",4);
  clean2Char(tag);
  printf("clean2Char test3: (%s)\n",tag);
  strncpy(tag,"CC",4);
  clean2Char(tag);
  printf("clean2Char test4: (%s)\n",tag);
  return 0;
}
 */
