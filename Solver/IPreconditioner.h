#pragma once

namespace ses {
	class IPreconditioner
	{
	public:
		IPreconditioner()
		{

		}
	};

	class DummyPreconditioner : public IPreconditioner {
	public:
		DummyPreconditioner()
		{

		}
	};
}
