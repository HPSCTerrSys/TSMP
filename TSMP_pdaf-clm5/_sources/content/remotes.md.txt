# TSMP remote repositories #

Relationship between the two remotes of the `TSMP` repository:

- Github
  ([https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP))
- Gitlab
  ([https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp](https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp))

Technically the two servers are two remote servers of the same
git-repository `TSMP`.

This means that both repositories have the same history of commits and
if you clone either remote repository, all branches that are part of
both repositories are identical. The difference is that the
Gitlab-remote contains additional branches, where TSMP-PDAF is
developed.

## Development scheme ##

- Github: `TSMP` general development
- Gitlab: `TSMP` is always kept up-to-date and changes are merged into
  the `TSMP-PDAF`-related branches.  
  
See also [development best practices](./best_practices)

## Privacy ##

- Github: Completely open.
- Gitlab: Only visible to invited members.


## More in-depth Example ##

Output of `TSMP`-git-repository's recent history including tags for
local branches and branches from the two remote repositories `origin`
(Gitlab) and `github` (Github):

``` text
*   05099a1 - (4 weeks ago)Merge branch 'master' into TSMP_pdaf_github - Johannes Keller (HEAD -> TSMP_pdaf, origin/TSMP_pdaf, github/TSMP_pdaf, TSMP_pdaf_github)
|\
| * 423ec0c - (4 weeks ago).gitignore: Ignore directories /cosmo4_21*/ - Johannes Keller (origin/master, github/master, github/HEAD, master-github, master)
* | 786909f - (4 weeks ago)Merge branch 'master' into TSMP_pdaf_github - Johannes Keller
|\|
| * b7beccd - (4 weeks ago)upgrade to ParFlow3.7 with possibilty of heterogeneous job submission - Abouzar Ghasemi
| * 8803904 - (4 weeks ago)upgrade to ParFlow3.7 with possibilty of heterogeneous job submission - Abouzar Ghasemi (tag: v1.3.3)
* | 7a26c8b - (5 weeks ago)Merge branch 'master-github' into TSMP_pdaf - Johannes Keller (tag: v1.2.3PDAF, origin/TSMP_pdaf_v1.2.3, TSMP_pdaf_v1.2.3)
|\|
| * 3cc6477 - (5 weeks ago)Correcting machine file for JURECA - Abouzar Ghasemi
| * 10a4cb3 - (6 weeks ago)Porting TSMP on JUSUF and JURECA_DC - Abouzar Ghasemi (tag: v1.2.3)
| * 8d54475 - (7 weeks ago)Machine files for JUWELS compatible with Stage2020 (Intel) and Stage2019a (Gnu) - Abouzar Ghasemi
* | 507714f - (8 weeks ago)compiler specific changes in src.Gnu / src.Intel in Oasis intf - Johannes Keller
```

Note that the commit `05099a1` is the current commit of the branches
`TSMP_pdaf, origin/TSMP_pdaf, github/TSMP_pdaf, TSMP_pdaf_github`. In
the local repository these are: 
- The local branch `TSMP_pdaf`
- the branch at remote origin (Gitlab) `origin/TSMP_pdaf`
- the branch at remote github (Github) `github/TSMP_pdaf`
- and the local version of the Github remote `TSMP_pdaf_github`
(introduced for automatic pushing and pulling, should always be at the
same commit as `TSMP_pdaf`).

The commit `423ec0c` is the current commit of `origin/master,
github/master, master-github, master` with similar meaning as before
for the `TSMP_pdaf` branches.

On commit `10a4cb3`, you see the tag `v1.2.3` for the version of
`TSMP` before ParFlow3.7 was introduced. From this version the commit
`7a26c8b` was derived with branches `origin/TSMP_pdaf_v1.2.3,
TSMP_pdaf_v1.2.3`. These branches only exist on Gitlab!
