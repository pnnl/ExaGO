# How to connect GitHub to a PNNL GitLab (push) Mirror

CI runners

path to ci.yml

in pnnl_mirror.yaml 
- generate 
	GIT_PASSWORD (done in gitLAB under settings>access tokens - give it write_repository perms + developer) 
	PNNL_PIPELINE_TRIGGER (done in gitLAB under CI/CD pipeline trigger)
- paste as variable in gitHUB


then github token for pnnl to push back results
GITHUB_CURL_HEADER


personal github settings - fine-grained tokens
	resource owner = pnnl
	select repositories = ppnl/exago
	permissions > repository permissions > commit statuses (read and write)
	copy this token
	go to gitLAB > settings > cicd > variables
		add variable
		type = file
		do not protect/mask/expand
		key = GITHUB_CURL_HEADER
		Value = Authorization: token <token value> 
https://ecp-ci.gitlab.io/docs/guides/build-status-gitlab.html 

	again for commit token
	permissions > repository permissions > contents (read and write)
	go to gitLAB > settings > cicd > variables
		add variable
		type = variable
		do not protect
		key = SPACK_GIT_TOKEN
		paste in value field