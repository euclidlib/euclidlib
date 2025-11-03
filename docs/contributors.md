# 🤝 Contributors

A heartfelt thank you to everyone who has contributed to **euclidlib**! 💫
Your ideas, work, and feedback have shaped this project into what it is today.
We recognize all kinds of contributions with the help of the [@all-contributors](https://allcontributors.org) bot.

---

```{code-cell} ipython3
from pathlib import Path
from IPython.display import Markdown, display

# Adjust the path if needed (assuming docs/ is next to README.md)
readme = Path("../README.md").read_text(encoding="utf-8")

start_tag = "<!-- ALL-CONTRIBUTORS-LIST:START -->"
end_tag = "<!-- ALL-CONTRIBUTORS-LIST:END -->"

if start_tag in readme and end_tag in readme:
    contributors = readme.split(start_tag)[1].split(end_tag)[0].strip()
    display(Markdown(contributors))
else:
    print("Contributors section not found in README.md")
```

💡 _Want to join them?_
Check out our [contribution guidelines](https://github.com/euclidlib/euclidlib) to learn how you can help — from improving documentation to adding new features or fixing bugs.
