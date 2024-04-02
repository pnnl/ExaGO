# How to connect GitHub to a PNNL GitLab (push) Mirror

To run on our HPC clusters at PNNL while hosting our code base on GitHub, we utilize PNNL's CI/CD services on GitLab. 

With GitLab premium, integration between GitHub and GitLab is included. The tier of GitLab that we have does not and did not initially support this. While switching to premium would alleviate some of the burden here, our technical requirements expanded having to test on more than one cluster. This is how we have architected a solution with the base level GitLab offering. 

The PNNL GitLab repository is a push mirror of the GitHub. So whenever a commit is pushed to a pull request or the `main` branch - the changes are pushed to the GitLab and a CI pipeline is triggered.

## Steps in linking the GitLab and GitHub

1. Push mirror & pipeline trigger

In `pnnl_mirror.yaml`, we use the variables `GIT_USER`, `GIT_PASSWORD`, and `PNNL_PIPELINE_TRIGGER`. `GIT_PASSWORD` and `PNNL_PIPELINE_TRIGGER` are generated in GitLab and then added to GitHub.

`GIT_USER` is the username of whoever will be authenticated when pushing to GitLab for the mirror action. Since we then manually trigger CI after that, and explicitly skip CI here, this name is cosmetic.

a) GIT_PASSWORD
Go to Settings > Access Tokens. Click `Add new token`.

Pick a reasonable name & expiration date. (ie "GITHUB_PUSH_PASSWORD")

Pick the `Developer` role.

Under `Select scopes`, select `write_repository`.

Create.

Go to GitHub > Settings > Secrets and variables > Actions. Click `New repository secret`.
Name it `GIT_PASSWORD` and paste in the value generated from GitLab.

Add secret.

b) PNNL_PIPELINE_TRIGGER
Go to Settings > CI/CD > Pipeline trigger tokens. Click `Add new token`.

Give it a name and click create. 

Go to GitHub > Settings > Secrets and variables > Actions. Click `New repository secret`.
Name it `PNNL_PIPELINE_TRIGGER` and paste in the value generated from GitLab.

2. Push back status & modules

In our module rebuild pipelines, we rebuild spack modules on each PNNL platform, then commit the new module paths/hashes back to the repository. 

a) Generate token #1 (`GITLAB_MIRROR_STATUS`)

Go to your GitHub profile > Settings > Developer Settings > Personal access tokens > Fine-grained tokens.

Generate new token 
	resource owner = pnnl
	select repositories = pnnl/exago
	permissions > repository permissions > commit statuses (read and write)
	copy this token

Go to GitLab > Settings > CI/CD > Variables
	add variable
	type = file
	do not protect/mask/expand
	key = GITHUB_CURL_HEADER
	Value = `Authorization: token <token value>`

See https://ecp-ci.gitlab.io/docs/guides/build-status-gitlab.html for more details.

b) Generate token #2 (`GITLAB_COMMIT`)

Go to your GitHub profile > Settings > Developer Settings > Personal access tokens > Fine-grained tokens.

Generate new token 
	resource owner = pnnl
	select repositories = pnnl/exago
	permissions > repository permissions > contents (read and write)
	copy this token

Go to GitLab > Settings > CI/CD > Variables
		add variable
		type = variable
		do not protect
		key = SPACK_GIT_TOKEN
		paste in value field


## Change path to `.gitlab-ci.yml`

Go to > Settings > CI/CD > General Pipelines, change the `CI/CD configuration file` to the correct path to `.gitlab-ci.yml`.

In our repo, the path is `.gitlab/pnnl/.gitlab-ci.yml`.
